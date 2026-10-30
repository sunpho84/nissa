#ifndef _NICE_SETUP_HPP
#define _NICE_SETUP_HPP

#include <cstddef>
#include <string>
#include <vector>

struct Evaluator;

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
  
  struct Correlations;
  
  /// Set up correlation builders and the legacy-input diagnostic function
  /// The specification must outlive the evaluator's registered functions
  void setupNice(Evaluator& evaluator,
                 Correlations& correlations);
}

#endif
