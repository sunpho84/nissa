#ifndef _NICE_SCRIPT_ANALYSIS_HPP
#define _NICE_SCRIPT_ANALYSIS_HPP

#include "nASTy.hpp"
#include "scriptCompletion.hpp"

#include <unordered_map>

/// AST and source ranges retained for a single editor revision
struct ParsedScript
{
  /// Owned original program pieces, including their source map
  ScriptSource source;
  
  /// Shared grammar tables also identify the expected terminals in diagnostics
  std::shared_ptr<const pp::Grammar> grammar;
  
  /// Complete AST or surviving AST fragments after bounded recovery
  std::vector<std::shared_ptr<ASTNode>> roots;
  
  /// Joined-source range for each node constructed by the normal AST actions
  std::unordered_map<const ASTNode*,SourceRange> ranges;
  
  /// Lexer/parser diagnostics, independent of Emacs presentation
  std::vector<ParseDiagnostic> diagnostics;
  
  /// Replacement extent when the requested position can complete an identifier
  std::optional<SourceRange> completion;
  
  /// Keywords obtained from the actual grammar's terminal definitions
  std::vector<std::string> keywords;
  
  /// Original lexer spans; no string views borrow a moved source string
  struct Token
  {
    /// Half-open original token extent
    SourceRange range;
    
    /// Whether this is a grammar identifier
    bool identifier;
    
    /// Whether this is the end-of-input sentinel
    bool end;
  };
  
  /// Lexer output retained for cursor lookup without tokenizing again
  std::vector<Token> tokens;
  
  /// Skipped text spans used to suppress completion inside comments
  std::vector<SourceRange> trivia;
  
  /// Invalid lexer spans, kept separate from recoverable parser diagnostics
  std::vector<SourceRange> lexicalErrors;
};

/// Marker spellings occur only in editor tokens, never in user-visible bindings
inline constexpr std::string_view completionMarker="\2";
inline constexpr std::string_view recoveryMarker="\1";

/// Parse with the actual lexer, LR tables, and AST actions, without evaluating script code
/// Recovery retains completed trees and inserts missing grammar terminals when possible
ParsedScript parseScriptForEditor(const ScriptSource& source,
                                  const std::optional<size_t> completionOffset={});

/// Resolve a cursor using retained lexer spans, without parsing or executing source
std::optional<SourceRange> scriptCompletionRange(const ParsedScript& script,
                                                const size_t byteOffset);

/// One editor revision shared by syntax checking and cursor-dependent AST visitors
struct ScriptAnalysisSession
{
  /// Canonical recovered AST, never marked with a completion cursor
  std::optional<ParsedScript> parsed;
  
  /// Bounded fallback for syntax which recovery could not retain at the cursor
  std::optional<ParsedScript> fallback;
  
  /// Source revision, advanced only when text or original-piece mapping changes
  size_t revision{0};
  
  /// Actual parse count, exposed for testing cache reuse
  size_t nParses{0};
  
  /// Replace a changed document; equal bytes with different source labels also invalidate
  void update(const ScriptSource& source)
  {
    bool same=parsed and parsed->source.text==source.text and parsed->source.pieces.size()==source.pieces.size();
    for(size_t i=0;same and i<source.pieces.size();i++)
      {
        const auto& old=parsed->source.pieces[i];
        const auto& next=source.pieces[i];
        same=old.begin==next.begin and old.size==next.size and old.origin==next.origin;
      }
    if(not same)
      {
        fallback.reset();
        parsed=parseScriptForEditor(source);
        revision++;
        nParses++;
      }
  }
  
  /// Query the retained AST; only lost recovery targets need a marked fallback parse
  ScriptCompletion complete(const Environment& environment,
                            const size_t byteOffset)
  {
    ScriptCompletion result;
    if(parsed)
      {
        result=completeScriptAt(environment,*parsed,byteOffset);
        const auto range=scriptCompletionRange(*parsed,byteOffset);
        if(not result.applicable and range)
          {
            if(not fallback or not fallback->completion or
               fallback->completion->begin!=range->begin or fallback->completion->end!=range->end)
              {
                fallback=parseScriptForEditor(parsed->source,byteOffset);
                nParses++;
              }
            result=completeScriptAt(environment,*fallback,byteOffset);
          }
      }
    
    return result;
  }
};

/// Write an escaped Lisp string as data, never executable code
void writeScriptLispString(std::ostream& output,
                           const std::string_view value);

/// Serialize all mapped diagnostics as (SOURCE BEGIN END LINE COLUMN MESSAGE) entries
std::string formatScriptDiagnosticsForEmacs(const ParsedScript& script);

#endif
