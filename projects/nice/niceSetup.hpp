#ifndef _NICE_SETUP_HPP
#define _NICE_SETUP_HPP

#include <cstddef>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

struct Evaluator;
struct ASTNodesNode;
struct ScriptCompletion;
struct ScriptSource;

namespace nissa::nice
{
  /// Query registered NICE metadata and preceding assignments without evaluating input or initializing NISSA
  ScriptCompletion completeNiceScript(const std::string_view text,
                                      const size_t byteOffset);
  
  /// Query assembled script pieces with original source locations and joined byte offsets
  ScriptCompletion completeNiceScript(const ScriptSource& source,
                                      const size_t byteOffset);
  
  /// Handle --complete-at OFFSET, --check-syntax, or --editor-server before MPI startup
  /// The server retains revisions and serves length-framed source pieces until input EOF
  /// Source text comes from input; offsets are local to it and output contains data only
  /// Ordered --context-file/--context-script options supply preceding pieces without execution
  /// --completion-source names the input piece; file errors are reported before MPI startup
  /// Syntax checks build a recovered AST and accept --check-source and --check-format text|emacs
  /// All diagnostics retain mapped source ranges; errors set exit status one without evaluating the AST
  /// Return no status for a normal launch, otherwise the query's exit status
  std::optional<int> runScriptQuery(const int narg,
                                        const char* const* arg,
                                        std::istream& input,
                                        std::ostream& output,
                                        std::ostream& diagnostic);
  
  /// One program source described on the command line
  struct ProgramPart
  {
    /// Distinguish directly supplied program text from a file path
    enum class Kind{Immediate,File};
    
    /// Interpretation of value
    Kind kind;
    
    /// Owned program text or the path supplied after -f
    std::string value;
    
    /// Position of the text or path operand in the original arg
    std::size_t argIndex;
  };
  
  /// Program sources in their original command-line order
  using LaunchArgs=
    std::vector<ProgramPart>;
  
  /// Describe program sources without reading files; MPI and rank information must be initialized
  /// Invalid launch options are reported through CRASH
  LaunchArgs parseLaunchArgs(const int narg,
                             const char* const* arg);
  
  /// Owned source text and its label for command-line or file diagnostics
  struct ProgramPiece
  {
    /// Source text retained while parsing and evaluating the program
    std::string text;
    
    /// Original argument index and, for file pieces, the supplied path
    std::string origin;
    
    /// Actual file path for editor navigation, empty for immediate arguments
    std::string filePath{};
  };
  
  /// Read files only on the master rank and broadcast ordered text and source labels
  /// All ranks must call this after MPI and rank initialization
  std::vector<ProgramPiece> readAndBroadcastProgramPieces(const LaunchArgs& parts);
  
  /// Join pieces with a terminating newline where needed, then parse one program
  /// Explicit blocks retain their scopes; documentation can cross piece boundaries
  /// The returned AST owns its text; retain it throughout evaluation
  ASTNodesNode parseProgramPieces(const std::vector<ProgramPiece>& pieces);
  
  struct PreparedProgram;
  
  /// Collectively load and evaluate a program, compile correlations, and retain deferred gauge recipes
  /// MPI and rank information must be initialized; the result owns all returned data
  PreparedProgram prepareProgram(const LaunchArgs& parts);
  
  struct RunOptions;
  
  struct RngLayout;
  struct CompiledCorrelations;
  struct GaugeRecipe;
  
  /// Collect each operation's double cost and validate volume-scaled scope budgets before execution
  /// Keep options alive and unchanged while using the returned layout
  RngLayout prepareRngLayout(const RunOptions& options,
                             const CompiledCorrelations& correlations,
                             const std::vector<GaugeRecipe>& gaugeRecipes);
  
  /// Register run-wide setters without performing geometry or numerical setup
  /// Options must outlive the registered functions; validate after program evaluation
  void registerRunOptions(Evaluator& evaluator,
                          RunOptions& options);
  
  struct Correlations;
  
  /// Set up correlation builders and the legacy-input diagnostic function
  /// The specification must outlive the evaluator's registered functions
  void setupNice(Evaluator& evaluator,
                 Correlations& correlations);
}

#endif
