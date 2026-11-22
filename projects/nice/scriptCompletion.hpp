#ifndef _NICE_SCRIPT_COMPLETION_HPP
#define _NICE_SCRIPT_COMPLETION_HPP

#include "nASTy.hpp"

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

/// Query a supplied scope and reconstruct completed assignments without executing input
/// Supports names, keywords, named arguments, aliases, and statically known expression members
/// Native calls and registered arithmetic use result metadata without invoking their callables
/// Brace scopes retain local visibility; unsupported expressions have unknown result types
/// Function parameters hide enclosing bindings but their types are not inferred
ScriptCompletion completeScriptAt(const Environment& environment,
                                  const std::string_view text,
                                  const size_t byteOffset);

/// Query assembled source while retaining original file/argument locations
/// The query offset and returned replacement bounds refer to joined-source bytes
ScriptCompletion completeScriptAt(const Environment& environment,
                                  const ScriptSource& source,
                                  const size_t byteOffset);

/// Serialize to Lisp data for an Emacs completion-at-point adapter; never executable Lisp
/// Shape: (BEGIN END ((NAME KIND DOCUMENTATION SOURCE LINE COLUMN OFFSET) ...)), or nil
std::string formatScriptCompletionForEmacs(const ScriptCompletion& completion);

#endif
