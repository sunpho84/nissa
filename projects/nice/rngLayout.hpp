#ifndef _NICE_RNG_LAYOUT_HPP
#define _NICE_RNG_LAYOUT_HPP

#include "runOptions.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <optional>
#include <vector>

namespace nissa::nice
{
  /// Consecutive doubles assigned to one operation within its execution scope
  struct RngAllocation
  {
    /// Execute once per gauge configuration or once per hit
    RngScope scope{RngScope::GaugeConf};
    
    /// First double within the scope's reserved range
    RngPosition offset{};
    
    /// Uniform-double consumption cost reported for one execution
    uint64_t doubleDrawCost{};
  };
  
  /// Prepared operation ranges sharing one volume-scaled budget per execution scope
  struct RngLayout
  {
    /// Borrowed run settings, kept alive and unchanged while using this layout
    const RunOptions& options;
    
    /// Consumption costs for gauge configuration preparation and one hit, accumulated at full width
    std::array<RngPosition,2> doubleDrawCosts{};
    
    /// Allocation for the optional origin at the beginning of each hit
    RngAllocation origin;
    
    /// Allocations in gauge recipe order, including zero-cost recipes
    std::vector<RngAllocation> gaugeOperations;
    
    /// Refer to the run settings before collecting any operation costs
    explicit RngLayout(const RunOptions& options) :
      options(options)
    {
      validate();
    }
    
    /// Number of doubles available to gauge configuration preparation or to any single hit
    RngPosition reservedDoubles() const
    {
      return RngPosition{{{options.physics.volume(),0,0,0}}}*options.physics.reservedRngDoublesPerSite;
    }
    
    /// Check a full-width total against the shared budget
    bool fits(const RngPosition& total) const
    {
      const RngPosition reserved=
        reservedDoubles();
      
      return not std::lexicographical_compare(std::rbegin(reserved.data),std::rend(reserved.data),
                                              std::rbegin(total.data),std::rend(total.data));
    }
    
    /// Append an operation's own double count to the selected scope
    RngAllocation addOperation(const RngScope scope,
                               const uint64_t doubleDrawCost)
    {
      auto& total=doubleDrawCosts[rngScopeIndex(scope)];
      const RngAllocation allocation{scope,total,doubleDrawCost};
      total+=doubleDrawCost;
      validate();
      
      return allocation;
    }
    
    /// Reject invalid geometry or costs exceeding either scope's reserved budget
    void validate() const
    {
      options.validate();
      
      for(size_t scope=0;scope<doubleDrawCosts.size();scope++)
        if(not fits(doubleDrawCosts[scope]))
          CRASH("NICE RNG %s operations exceed reserved doubles per site",
                scope==0?"gauge configuration":"hit");
    }
    
    /// Position an operation without advancing earlier gauge configurations, hits, or operations
    /// An absent hit selects gauge configuration preparation; a present hit selects slot iHit+1
    /// The selected execution scope must match the operation allocation
    /// A separate 64-bit gauge configuration index keeps positions independent of requested hit counts
    OperationRngRange operationRange(const Sitmo::Rng& rng,
				     const RngAllocation& operation,
				     const uint64_t iGaugeConf,
				     const std::optional<uint64_t> iHit=std::nullopt) const
    {
      rngScopeIndex(operation.scope);
      const RngScope scope=iHit?RngScope::Hit:RngScope::GaugeConf;
      
      if(operation.scope!=scope)
        CRASH("NICE RNG hit index does not match the operation's scope");
      
      if(not fits(operation.offset+operation.doubleDrawCost))
        CRASH("NICE RNG operation exceeds its scope's reserved range");
      
      const uint64_t slot=
        iHit?checkedAdd(*iHit,uint64_t{1},"NICE RNG hit index overflow"):0;
      const RngPosition position=
        RngPosition{{{slot,iGaugeConf,0,0}}}*options.physics.volume()*options.physics.reservedRngDoublesPerSite+operation.offset;
      
      return {rng,position,operation.doubleDrawCost};
    }
  };
}

#endif
