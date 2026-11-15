#ifndef _NICE_WALL_TIME_HPP
#define _NICE_WALL_TIME_HPP

#include "routines/ios.hpp"
#include "runOptions.hpp"
#include "routines/mpi_routines.hpp"

#include <charconv>
#include <chrono>
#include <cstdlib>
#include <string>
#include <string_view>

namespace nissa::nice
{
  /// Collective admission of new work using IB's pessimistic completed-hit average
  struct WallTimeBudget
  {
    /// Monotonic clock used to measure elapsed execution time
    using Clock=std::chrono::steady_clock;
    
    /// Slurm variable containing the projected allocation end as a Unix timestamp
    static constexpr char slurmEndTimeVariable[]{"SLURM_JOB_END_TIME"};
    
    /// Start of this invocation, including library initialization and script preparation
    Clock::time_point start;
    
    /// Allowed seconds since start, resolved once on the master rank
    double limitSeconds{};
    
    /// Fixed allowance added to the hit estimate
    double safetySeconds;
    
    /// Hits computed and committed in this invocation, excluding restored checkpoints
    size_t completedHits{0};
    
    /// Resolve explicit duration, Slurm end timestamp, or the one-minute fallback
    /// Supplied clock readings allow deterministic validation without sleeping
    WallTimeBudget(const RunOptions& options,
                   const Clock::time_point start,
                   const char* slurmEnd=std::getenv(slurmEndTimeVariable),
                   const Clock::time_point now=Clock::now(),
                   const double unixNow=std::chrono::duration<double>(
                     std::chrono::system_clock::now().time_since_epoch()).count()) :
      start(start),
      safetySeconds(options.wallTimeSafetySeconds)
    {
      std::string source;
      
      if(is_master_rank())
        {
          if(options.wallTimeSeconds)
            {
              limitSeconds=*options.wallTimeSeconds;
              source="script";
            }
          else if(slurmEnd!=nullptr and *slurmEnd!='\0')
            {
              const std::string_view text(slurmEnd);
              uint64_t timestamp{};
              const auto parsed=std::from_chars(text.data(),text.data()+text.size(),timestamp);
              
              if(parsed.ec!=std::errc{} or parsed.ptr!=text.data()+text.size() or timestamp==0)
                CRASH("NICE invalid %s '%s'",slurmEndTimeVariable,slurmEnd);
              
              limitSeconds=elapsed(now)+static_cast<double>(timestamp)-unixNow;
              source=slurmEndTimeVariable;
            }
          else
            {
              limitSeconds=60.0;
              source="one-minute fallback";
            }
        }
      
      limitSeconds=broadcast(limitSeconds);
      MASTER_PRINTF("NICE wall time: %s, limit %.17g s from startup, safety %.17g s\n",
                    source.c_str(),limitSeconds,safetySeconds);
    }
    
    /// Seconds since startup, including preparation, checkpoint I/O, and finalization
    double elapsed(const Clock::time_point now=Clock::now()) const
    {
      return std::chrono::duration<double>(now-start).count();
    }
    
    /// Require room for the IB average plus ten percent and the explicit allowance
    /// Before the first completed hit only the deadline and allowance are known
    bool allowsWork(const Clock::time_point now=Clock::now()) const
    {
      bool allowed=false;
      
      if(is_master_rank())
        {
          const double passed=elapsed(now);
          const double remaining=limitSeconds-passed;
          const double estimate=completedHits==0?0.0:1.1*passed/completedHits;
          allowed=remaining>estimate+safetySeconds;
          MASTER_PRINTF("NICE wall time: remaining %.17g s, hit estimate %.17g s (%zu completed), safety %.17g s: %s\n",
                        remaining,estimate,completedHits,safetySeconds,
                        allowed?"continue":"stop");
        }
      
      return broadcast(allowed);
    }
    
    /// Add a hit only after its checkpoint has been committed
    void recordCompletedHit()
    {
      completedHits++;
    }
  };
}

#endif
