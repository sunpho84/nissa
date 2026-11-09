#ifndef _NICE_RNG_HPP
#define _NICE_RNG_HPP

#include "base/randomDistributions.hpp"
#include "base/sitmo.hpp"
#include "geometry/geometry_lx.hpp"

#include <cstdint>

namespace nissa::nice
{
  /// Execution scope selecting gauge configuration-wide or per-hit draws
  enum class RngScope{GaugeConf,Hit};
  
  /// Return a checked index for the two execution scopes
  inline size_t rngScopeIndex(const RngScope scope)
  {
    const size_t index=static_cast<size_t>(scope);
    
    if(index>=2)
      CRASH("NICE RNG operation has an invalid scope");
    
    return index;
  }
  
  /// Four-word position in Sitmo's indexed random stream
  using RngPosition=
    Sitmo::Counter;
  
  /// Bounded random access to the doubles assigned to one operation execution
  struct OperationRngRange
  {
    /// Generator snapshot positioned at this operation's first integer draw
    Sitmo::Rng rng;
    
    /// Beginning of this operation's allocation in double units
    RngPosition position;
    
    /// Number of accessible doubles reported by the operation
    uint64_t nDoubles{};
    
    /// Position the copied generator once, preserving the parent's state
    OperationRngRange(const Sitmo::Rng& runRng,
		   const RngPosition& position,
		   const uint64_t nDoubles) :
      rng(runRng),
      position(position),
      nDoubles(nDoubles)
    {
      rng.state+=position*uint64_t{2};
    }
    
    /// Convert a local double index to its 65-bit offset in integer draws
    CUDA_HOST_AND_DEVICE
    RngPosition integerOffset(const uint64_t iDouble) const
    {
      if(iDouble>=nDoubles)
        CRASH("NICE RNG draw exceeds the operation's reported double count");
      
      return RngPosition{{{iDouble<<1,iDouble>>63,0,0}}};
    }
    
    /// Generate one uniform double independently of execution order or MPI layout
    CUDA_HOST_AND_DEVICE
    double uniformAt(const uint64_t iDouble) const
    {
      RngView view(rng,integerOffset(iDouble));
      
      return uniformDouble(view);
    }
  };
  
  /// Implicit per-hit operation selecting the shared source origin
  struct GlobalOriginOperation
  {
    /// Generate a random origin when enabled, otherwise use the lattice origin
    bool enabled{true};
    
    /// Report the double count for one hit
    uint64_t rngDoubles() const
    {
      return enabled?uint64_t{NDIM}:uint64_t{0};
    }
    
    /// Generate the origin from the assigned doubles without scaling its allocation by the volume
    Coords execute(const OperationRngRange& range) const
    {
      if(range.nDoubles!=rngDoubles())
        CRASH("NICE origin operation does not match its RNG allocation");
      
      Coords origin{};
      
      if(enabled)
        for(int mu=0;mu<NDIM;mu++)
          origin[mu]=static_cast<int>(range.uniformAt(mu)*glbSize[mu]);
      
      return origin;
    }
  };
}

#endif
