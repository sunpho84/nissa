#ifndef _NICE_RUN_POLICY_HPP
#define _NICE_RUN_POLICY_HPP

#include "base/debug.hpp"

#include <cmath>
#include <cstdint>
#include <optional>
#include <string>

namespace nissa::nice
{
  /// Execution and output policy, independent of the identity of one numerical hit
  struct RunPolicy
  {
    /// Number of stochastic hits per gauge configuration
    int nHits{1};
    
    /// Retain completed hit files and compatibility metadata for later hit-count increases
    bool keepCheckpoints{true};
    
    /// Filename suffix used to keep independent runs in one output directory
    std::string outputSuffix;
    
    /// Optional stop-request file, interpreted relative to the process working directory
    std::string stopFile;
    
    /// Seconds between running-file updates; a file older than twice this can be reclaimed
    int runningFileUpdateSeconds{60};
    
    /// Explicit seconds from program startup, otherwise use Slurm or a one-minute fallback
    std::optional<double> wallTimeSeconds;
    
    /// Extra seconds reserved beyond the pessimistic completed-hit estimate
    double wallTimeSafetySeconds{0.0};
    
    /// Maximum consecutive attempts without a committed hit; zero permits unlimited retries
    int maxTrials{3};
    
    /// Require meaningful work limits, timing values, and output names
    void validate() const
    {
      if(nHits<=0)
        CRASH("NICE run requires positive nHits");
      
      if(maxTrials<0)
        CRASH("NICE maximum trials must be nonnegative");
      
      if(runningFileUpdateSeconds<=0)
        CRASH("NICE running-file update interval must be positive");
      
      if(wallTimeSeconds and (not std::isfinite(*wallTimeSeconds) or *wallTimeSeconds<=0))
        CRASH("NICE wall time must be finite and positive");
      
      if(not std::isfinite(wallTimeSafetySeconds) or wallTimeSafetySeconds<0)
        CRASH("NICE wall-time safety allowance must be finite and nonnegative");
      
      for(const unsigned char character : outputSuffix)
        if(not ((character>='a' and character<='z') or
                (character>='A' and character<='Z') or
                (character>='0' and character<='9') or
                character=='_' or character=='-'))
          CRASH("NICE output suffix may contain only letters, digits, '_' and '-'");
    }
  };
}

#endif
