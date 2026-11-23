#ifndef _NICE_SCRIPT_SESSION_HPP
#define _NICE_SCRIPT_SESSION_HPP

#include <iosfwd>

struct Environment;

/// Serve length-framed editor queries against a persistent AST without executing scripts
/// Each request supplies ordered labelled source pieces; equal revisions reuse their AST
/// Responses are one-line Lisp data: (ID STATUS REVISION NPARSES PAYLOAD MESSAGE)
int runScriptSession(const Environment& environment,
                     std::istream& input,
                     std::ostream& output,
                     std::ostream& diagnostic);

#endif
