#ifndef _NICE_RUN_OPTIONS_HPP
#define _NICE_RUN_OPTIONS_HPP

#include "base/debug.hpp"
#include "rng.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/math_routines.hpp"

#include <cstdint>
#include <limits>

namespace nissa::nice
{
  /// Run-wide settings recorded before any geometry or numerical setup
  struct RunOptions
  {
    /// Temporal fermion boundary, independent of spatial momentum
    enum class TemporalBoundary{Periodic,Antiperiodic};
    
    /// Common spatial extent, required before execution
    int L{};
    
    /// Temporal extent, required before execution
    int T{};
    
    /// Seed inherited from Sitmo unless changed by setRun or a C++ caller
    int64_t seed{Sitmo::Rng::DEFAULT_SEED()};
    
    /// Number of stochastic hits per gauge configuration
    int nHits{1};
    
    /// Implicit origin-selection operation performed once per hit
    GlobalOriginOperation drawGlobalOriginTask;
    
    /// Independent copies per hit, kept separate when averaging hits
    int nCopies{1};
    
    /// Source spin dilution, disabled by default
    bool diluteSpin{false};
    
    /// Source color dilution, disabled by default
    bool diluteColor{false};
    
    /// Fermion temporal boundary, normally antiperiodic
    TemporalBoundary temporalBoundary{TemporalBoundary::Antiperiodic};
    
    /// Uniform doubles reserved per site for gauge configuration preparation and separately for each hit
    uint64_t reservedRngDoublesPerSite{1'000'000};
    
    /// Compute the global volume, rejecting invalid geometry before arithmetic
    uint64_t volume() const
    {
      if(L<=0 or T<=0)
        CRASH("NICE geometry requires positive L and T");
      
      uint64_t result=static_cast<uint64_t>(T);
      
      for(int mu=1;mu<NDIM;mu++)
        result=checkedMultiply(result,static_cast<uint64_t>(L),"NICE geometry volume overflow");
      
      // initGrid accumulates volume in int64_t
      if(result>static_cast<uint64_t>(std::numeric_limits<int64_t>::max()))
        CRASH("NICE geometry volume exceeds the signed 64-bit range");
      
      return result;
    }
    
    /// Validate the effective settings after all program pieces have executed
    void validate() const
    {
      volume();
      
      if(reservedRngDoublesPerSite==0)
        CRASH("NICE RNG requires a positive number of reserved doubles per site");
      
      if(seed<0 or seed>std::numeric_limits<uint32_t>::max())
        CRASH("NICE seed must fit Sitmo uint32_t range");
      
      if(nHits<=0 or nCopies<=0)
        CRASH("NICE run requires positive nHits and nCopies");
      
      if(temporalBoundary!=TemporalBoundary::Periodic and
         temporalBoundary!=TemporalBoundary::Antiperiodic)
        CRASH("NICE run has an invalid temporal boundary");
    }
  };
}

#endif
