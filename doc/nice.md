# NICE correlation preparation

NICE provides C++ types and nASTy bindings for describing correlations,
turning them into validated tasks, ordering those tasks, and writing a legacy
prepareInput representation. The correlation types live in
[`projects/nice/`](../projects/nice/). They do not depend on the old
prepareInput adapter. Program loading, run-setting bindings, and configuration
list loading are available as library calls;
a standalone frontend and numerical execution remain to be implemented.

NICE implementation work follows the
[`implementation conventions`](nissa-conventions.md).

## Current data flow

```text
C++ construction or nASTy script
              |
              v
correlation parts and Correlations
              |
              v
buildTaskGraph -> CorrelationTaskGraph
              |
              v
scheduleCorrelations -> CompiledCorrelations
              |
              v
writeLegacyInput -> prepareInput-compatible text
```

The types through `CompiledCorrelations` are independent of the legacy text
format. `writeLegacyInput` is a diagnostic and migration output. It does not
execute correlations.

## Program input and run settings

[`niceSetup.hpp`](../projects/nice/niceSetup.hpp) provides command-line parsing,
master-rank program-file reading and broadcast, and top-level AST assembly.
`parseLaunchArgs(narg,arg)` preserves positional statements and `-f PATH`
pieces in order; `--` ends option parsing. Call it after MPI and rank setup.
`readAndBroadcastProgramPieces()` returns owned source strings and labels.
`parseProgramPieces()` parses them separately and merges only their top-level
statement lists. Explicit braces retain an inner scope. Keep the source
pieces unchanged and alive alongside the AST throughout evaluation.
`prepareProgram(parts)` connects loading, evaluation, run-setting validation,
correlation compilation, and configuration list loading. Its
`PreparedProgram` result owns `options`, `correlations`, and
`configurationEntries`; it does not borrow the local AST, evaluator, or source text.
After MPI and rank setup, callers including `nice.hpp` and `niceSetup.hpp` can
use:

```cpp
const auto prepared=prepareProgram(parseLaunchArgs(narg,arg));
const auto legacyText=writeLegacyInput(prepared.correlations);
```

This prepares instructions without initializing the grid or running numerical
operations. Source-aware syntax diagnostics and a standalone executable remain pending.

[`RunOptions`](../projects/nice/runOptions.hpp) records caller-owned run-wide
settings. Call `registerRunOptions(evaluator,options)` alongside
`setupNice(evaluator,correlations)` to expose:

```text
setGeometry(.L=24,.T=48);
setRun(.seed=1234,.nHits=8,.nCopies=2);
setDilution(.diluteSpin=1,.diluteColor=1);
setTemporalBoundary("antiperiodic");
```

Each setter takes the complete argument group shown above. Later calls
replace that group's previous settings, including calls in later command-line
pieces. They only record settings: they do not initialize geometry, seed the
RNG, allocate fields, or perform numerical work. The `RunOptions` object must
outlive its registered functions. Call `options.validate()` after evaluation
and before setup. C++ callers can populate the same aggregate directly.

Both geometry extents, hit count, and copy count must be positive. The total
volume must fit the signed 64-bit representation used by NISSA geometry setup;
this catches accidental oversized extents before `initGrid`. The seed defaults
to Sitmo's `3472291050` and must fit its `uint32_t` range. Hits and copies
default to one. Spin and color dilution default to disabled; `setDilution`
accepts only `0` or `1`. The temporal boundary defaults to antiperiodic and
also accepts `"periodic"`. `setRun` currently requires all three arguments;
nASTy integer literals and its seed argument are `int`, so script overrides
above `INT_MAX` need future support. Copies remain separate
from hits and must not be averaged as extra hits. Numerical interpretation of
these settings remains part of driver development.

## Configuration lists

`configurationsFromFile("confs.list")` registers a list during script
evaluation. Each nonblank, noncomment line contains a configuration path, one
or more spaces or one tab, and an output directory. When either path contains
spaces, use exactly one tab to separate the two complete paths. Without a tab,
the line must contain exactly two whitespace-separated paths. A line whose
first nonwhitespace character is `#` is a comment. There is no entry count.
Multiple calls append entries in script order.

The master rank reads each list after evaluation and broadcasts its path
and contents. Every rank parses the same ordered rows. Relative paths in a row
are resolved against the list directory; the list path itself is
resolved from the working directory. Duplicate output directories and malformed
rows are rejected before numerical work. Input configuration files are checked
only when the driver loads them. The resulting entries are owned by
`PreparedProgram::configurationEntries`. Gauge preparation and other
configuration-source forms remain future work.

## Correlation parts and requests

[`corrParts.hpp`](../projects/nice/corrParts.hpp) and
[`corrParts.cpp`](../projects/nice/corrParts.cpp) define `SpatMomentum`, `Smear`,
`Phase`, `Gamma`, `DeltaT`, `Prop`, `Source`, `Oper`, and `Line`. Products build
line expressions; line addition and scalar multiplication build weighted
combinations. The types can be used directly from C++ without nASTy.
`Gamma` supports only the 16 standard matrices, indexed 0–15, for line
extensions and contractions; index 0 is the identity.
`Prop` remains an aggregate: omitted members copy the corresponding values
from the process-wide `Prop::defaults` at construction. Initially `kappa` is
`0.125` and the other values are zero. The Wilson sign `r` must be `-1` or
`+1` before graph validation or conjugation; zero means unset.
`setDefaultPropR(-1)` changes the sign inherited by later propagators;
`setDefaultPropR(0)` clears that default. Explicit constructor arguments
always win, and existing objects retain their values. Defaults are shared
by C++ callers and all evaluators in the process.
`setDefaultPropKappa(kappa)` and `setDefaultPropMass(mass)` change the
corresponding propagator defaults. `Smear` uses the same rules through
`Smear::defaults`, initially zero for kappa, iteration count, and momentum.
Use `setDefaultSmearKappa(kappa)` and `setDefaultSmearN(n)` to change its
strength and count. Kappa and mass setters take real arguments (for example
`0.0`); the iteration setter takes an integer argument (for example `20`).
Counts must be finite, nonnegative integers representable
by the backend's `int`; graph validation and the count setter reject other
values, and the legacy writer emits an integer count.

```cpp
setDefaultPropKappa(0.125);
setDefaultPropMass(0.01);
setDefaultPropR(-1);
setDefaultSmearKappa(0.3);
setDefaultSmearN(20);
prop=Prop();
smear=Smear();
unsmeared=Smear(.n=0);
```

`SpatMomentum` has three spatial components in both C++ and scripts. Temporal
behavior is specified separately through a boundary condition, normally
antiperiodic.

`dag()` denotes the Hermitian adjoint and reverses operator-product order.
`Prop::dag()` flips the Wilson sign multiplying the imaginary gamma-five term
and reverses momentum; `Phase::dag()` reverses its phase momentum. Standard
gamma indices 0–5 are Hermitian and 6–15 are anti-Hermitian. `Gamma::dag()`
and `dag(Pars)` therefore return an `Oper`, retaining the adjoint sign in its
real scalar weight. Existing real weights are unchanged by conjugation.

[`correlations.hpp`](../projects/nice/correlations.hpp) defines `Correlations`
as a collection of named `Correlation` blocks. Each block holds distinct
sink/source gamma pairs and ordered backward/forward line pairs. Block names
must be nonempty and unique. Duplicate trace pairs within a block are
rejected. A `Source` has a realization identity: copying it retains that
identity, while constructing another source allocates a new one. Scripts can
read its `tSelect` field but cannot assign its internal `count` identity.

A C++ caller can use [`nice.hpp`](../projects/nice/nice.hpp) to own requests
and compile them:

```cpp
using namespace nissa::nice;

Nice nice;
auto& block=nice.correlations.addBlock("pion");
block.addGammas(5,5);
const Line line(Source{0});
block.addTrace(line,line);

const CompiledCorrelations compiled=nice.compile();
const std::string text=writeLegacyInput(compiled);
```

## Lowering, validation, and scheduling

[`graphCompiler.hpp`](../projects/nice/graphCompiler.hpp) lowers nested line
expressions into a [`CorrelationTaskGraph`](../projects/nice/correlationTaskGraph.hpp).
Its tasks create sources, extend or select lines, and contract ordered line
pairs. Dependencies carry coefficients. Equivalent operations and inputs
share a task; user labels and block names do not determine task identity.
Tasks receive stable IDs in dependency-first order. Output blocks retain
request order and refer to their trace task IDs.

Graph construction validates task IDs, dependency order, input shapes,
coefficients, labels, and output references. Extension and contraction inputs
must come from sources or line extensions; contraction results cannot be used
as fields. Call `validate()` after directly
editing the public graph fields. A graph may be scheduled more than once.

[`schedulingStrategies.hpp`](../projects/nice/schedulingStrategies.hpp)
provides `scheduleByTaskId` and `scheduleByOutput`. The former follows the
graph's stable ID order. The latter visits each requested trace and its
dependencies before moving to the next trace, then includes any remaining
tasks. A strategy must return every graph task ID exactly once, including
unused tasks. `scheduleCorrelations` applies a strategy and returns
[`CompiledCorrelations`](../projects/nice/compiledCorrelations.hpp), an ordered
instruction list plus output references. Its validation requires every task
exactly once, with dependencies preceding users and outputs referencing valid
contractions. Scheduling changes order; it does not change the graph's
operations.

## Script bindings

[`niceSetup.cpp`](../projects/nice/niceSetup.cpp) registers the correlation
parts, their supported native binary operations (including `Line + Line`),
and these host functions:

- `setDefaultPropKappa(kappa)` and `setDefaultPropMass(mass)` set propagator defaults
- `setDefaultSmearKappa(kappa)` and `setDefaultSmearN(n)` set smearing defaults
- `setDefaultPropR(r)` sets the default Wilson sign (`-1` or `+1`), or clears it with `0`
- `corrBlock(name)` creates a named block
- `corrGammas(name,sink,source)` adds a gamma pair to that block
- `corrTrace(name,backward,forward)` adds an ordered line pair
- `writeLegacyInput()` compiles the current requests with the stable ID
  strategy and returns the legacy text

For example, with an `Evaluator` configured by `setupNice(evaluator,
correlations)`, a script can express:

```text
{
  corrBlock(.name="pion");
  corrGammas(.name="pion",.sink=5,.source=5);
  src=Source(0);
  backward=Prop(.kappa=0.12,.r=1)*src;
  forward=Gamma(1)*backward;
  corrTrace(.name="pion",.backward=backward,.forward=forward);
  print(writeLegacyInput());
}
```

The caller keeps the `Correlations` object, or its owning `Nice` instance,
alive while the evaluator's registered functions use it.
[`nASTy.hpp`](../projects/nice/nASTy.hpp)
provides generic host-object registration and positional/named argument
conversion; [`host-objects.md`](../projects/nice/host-objects.md) documents
those mechanisms.

The parser accepts multiple top-level statements without enclosing braces.
Braces create a nested block. `parseProgramPieces()` combines the
top-level statements from files and command-line text into one outer scope;
names declared inside an explicit block remain local to that block.

Binary operator precedence follows C/C++ for the supported operators.
Operators in each group share a precedence level: assignments, equality
comparisons, relational comparisons, addition/subtraction, and
multiplication/division/remainder. Arithmetic and comparisons associate to
the left; assignments associate to the right. Thus `8-3+2` is `(8-3)+2`,
`8*3/2` is `(8*3)/2`, and `a=b=c` is `a=(b=c)`.

Logical operators evaluate the left operand first. `and`/`&&` skip the right
operand when the left is false; `or`/`||` skip it when the left is true.
Their results are Boolean, using the same conversion rules as `if` and `for`.

Compound assignments (`+=`, `-=`, `*=`, `/=`) update script variables, list
elements, and struct members using native C++ assignment conversions.
Registered host arithmetic also supports the corresponding compound
assignments on script variables, for operand pairs supported by the concrete
C++ types and registered with the evaluator. Plain assignment may replace
a script variable's stored type.

Variadic script parameters collect extra arguments in `__VA__ARGS__`.
Use `...` to collect values or `&...` to collect references; each argument
in a reference pack must be passed explicitly by reference. Fixed parameters
may precede either kind of pack, and an empty pack is allowed:

```cpp
fun first(...){ return __VA__ARGS__[0]; }
fun setFirst(&...){ __VA__ARGS__[0]=7; }
x=1;
setFirst(&x);
```

## Legacy output and verification

[`legacyInputWriter.cpp`](../projects/nice/legacyInputWriter.cpp) renders the
compiled instructions as prepareInput source, propagator, and contraction
tables. It derives display names locally, including deterministic generated
names, and represents weighted dependencies with the legacy `LINCOMB` field.
Those names do not affect task identity or scheduling. The writer rejects
whitespace, the reserved name `LINCOMB` in any letter case, and names beginning
with `/*` in source labels, line labels, and output block names. A fixture in
[`validate.cpp`](../projects/nice/validate.cpp) compares the output byte for
byte with the old `Run` printer.

With ParsePact enabled, build and run the NICE validation executable using:

```sh
make -C <build-directory>/projects check
```

The single validation executable checks correlation parts, graph validation,
lowering, scheduling, legacy output, the script binding chain, the `Nice`
facade, and generic host-object registration. The old
`prepareInput.cpp` adapter and its `--test-prepare-*` modes have been removed.

