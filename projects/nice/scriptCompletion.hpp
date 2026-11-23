#ifndef _NICE_SCRIPT_COMPLETION_HPP
#define _NICE_SCRIPT_COMPLETION_HPP

#include "nASTy.hpp"

struct ParsedScript;

/// Documented replacement offered by the non-executing completion query
struct ScriptCompletionCandidate
{
  /// Text inserted in place of the identifier being completed
  std::string name;
  
  /// Entity, keyword, parameter, or member classification
  EntityInfo info;
  
  /// Formatted help, without evaluating defaults or calling script code
  std::string documentation;
};

/// Completion replacement range and candidates for one editor query
struct ScriptCompletion
{
  /// Whether identifier completion applies at this position
  bool applicable{false};
  
  /// Zero-based byte offset of the replacement start in the supplied text
  size_t begin{0};
  
  /// Exclusive zero-based byte offset of the replacement end
  size_t end{0};
  
  /// Prefix-matching candidates in alphabetical order
  std::vector<ScriptCompletionCandidate> candidates;
};

/// Query a supplied scope through the actual parser and a static AST visitor, without executing input
/// Supports names, keywords, named arguments, aliases, and statically known expression members
/// Native calls and registered arithmetic use result metadata without invoking their callables
/// Brace scopes retain local visibility; unsupported expressions have unknown result types
/// Calls specialize script parameters and struct members using argument/default AST metadata
/// Parameters inside an unspecialized declaration still hide enclosing bindings with unknown types
ScriptCompletion completeScriptAt(const Environment& environment,
                                  const std::string_view text,
                                  const size_t byteOffset);

/// Query assembled source while retaining original file/argument locations
/// The query offset and returned replacement bounds refer to joined-source bytes
ScriptCompletion completeScriptAt(const Environment& environment,
                                  const ScriptSource& source,
                                  const size_t byteOffset);

/// Visit a retained editor AST without reparsing its source
/// Ordinary snapshots support different cursor positions through retained lexer spans
/// Marked recovery snapshots support offsets within their marked identifier
ScriptCompletion completeScriptAt(const Environment& environment,
                                  const ParsedScript& parsed,
                                  const size_t byteOffset);

/// Serialize to Lisp data for an Emacs completion-at-point adapter; never executable Lisp
/// Shape: (BEGIN END ((NAME KIND DOCUMENTATION SOURCE LINE COLUMN OFFSET) ...)), or nil
std::string formatScriptCompletionForEmacs(const ScriptCompletion& completion);

/// Serialize NUL-separated bounds and shell replacements, adding a call for a unique callable
/// Required parameters of its first signature become placeholders; optional/variadic ones are omitted
/// Preserve existing argument lists and identifiers completed with point inside their spelling
std::string formatScriptCompletionForBash(const ScriptCompletion& completion,
                                         const std::string_view text,
                                         const size_t byteOffset);

#endif
