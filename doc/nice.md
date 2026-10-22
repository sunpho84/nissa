# NICE correlation preparation

NICE provides C++ types and nASTy bindings for describing correlations,
turning them into validated tasks, ordering those tasks, and writing a legacy
prepareInput representation. The correlation types live in
[`projects/nice/`](../projects/nice/). They do not depend on the old
prepareInput adapter. A program frontend and numerical execution are the
next development steps.

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

## Correlation parts and requests

[`corrParts.hpp`](../projects/nice/corrParts.hpp) and
[`corrParts.cpp`](../projects/nice/corrParts.cpp) define `SpatMomentum`, `Smear`,
`Phase`, `Gamma`, `DeltaT`, `Prop`, `Source`, `Oper`, and `Line`. Products build
line expressions; line addition and scalar multiplication build weighted
combinations. The types can be used directly from C++ without nASTy.
`SpatMomentum` has three spatial components in both C++ and scripts. Temporal
behavior is specified separately through a boundary condition, normally
antiperiodic.

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
coefficients, labels, and output references. Call `validate()` after directly
editing the public graph fields. A graph may be scheduled more than once.

[`schedulingStrategies.hpp`](../projects/nice/schedulingStrategies.hpp)
provides `scheduleByTaskId` and `scheduleByOutput`. The former follows the
graph's stable ID order. The latter visits each requested trace and its
dependencies before moving to the next trace, then includes any remaining
tasks. `scheduleCorrelations` applies a strategy and returns
[`CompiledCorrelations`](../projects/nice/compiledCorrelations.hpp), an ordered
instruction list plus output references. Its validation requires every task
exactly once, with dependencies preceding users and outputs referencing valid
contractions. Scheduling changes order; it does not change the graph's
operations.

## Script bindings

[`niceSetup.cpp`](../projects/nice/niceSetup.cpp) registers the correlation
parts, their supported native binary operations (including `Line + Line`),
and these host functions:

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
Braces create a nested block. A future program frontend can combine the
top-level statements from files and command-line text into one outer scope;
names declared inside an explicit block remain local to that block.

## Legacy output and verification

[`legacyInputWriter.cpp`](../projects/nice/legacyInputWriter.cpp) renders the
compiled instructions as prepareInput source, propagator, and contraction
tables. It derives display names locally, including deterministic generated
names, and represents weighted dependencies with the legacy `LINCOMB` field.
Those names do not affect task identity or scheduling. A fixture in
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

