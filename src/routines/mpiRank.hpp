#ifndef _MPITYPES_HPP
#define _MPITYPES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <expr/comp.hpp>

namespace nissa
{
  DECLARE_DYNAMIC_COMP(MpiRankCoord);
  
  DECLARE_DYNAMIC_COMP(MpiRank);
  
  /// Master rank
  constexpr MpiRank masterRank=0;
  
  PROVIDE_RESOURCE(thisRank,MpiRank);
  
  /// Check if this is the master rank
  INLINE_FUNCTION
  bool isMasterRank()
  {
    return thisRank==masterRank;
  }
}

#endif
