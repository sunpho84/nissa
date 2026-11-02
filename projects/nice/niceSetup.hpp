#ifndef _NICE_SETUP_HPP
#define _NICE_SETUP_HPP

#include <cstddef>
#include <string>
#include <vector>

struct Evaluator;
struct ASTNodesNode;

namespace nissa::nice
{
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
  };
  
  /// Read files only on the master rank and broadcast ordered text and source labels
  /// All ranks must call this after MPI and rank initialization
  std::vector<ProgramPiece> readAndBroadcastProgramPieces(const LaunchArgs& parts);
  
  /// Parse each piece and merge top-level statements, retaining explicit block scopes
  /// Keep pieces unchanged and alive, and retain the returned AST throughout evaluation
  ASTNodesNode parseProgramPieces(const std::vector<ProgramPiece>& pieces);
  
  struct PreparedProgram;
  
  /// Collectively load, evaluate, validate settings, and compile a NICE program
  /// MPI and rank information must be initialized; the result owns all returned data
  PreparedProgram prepareProgram(const LaunchArgs& parts);
  
  struct RunOptions;
  
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
