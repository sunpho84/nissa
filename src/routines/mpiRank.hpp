#ifndef _MPIRANK_HPP
#define _MPIRANK_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/universe.hpp>
#include <expr/comp.hpp>

namespace nissa
{
  DECLARE_DYNAMIC_COMP(MpiRankCoord);
  
  /// Mpi Rank coordinates
  using MpiRankCoords=
    Coords<MpiRankCoord>;
  
  namespace resources
  {
    inline MpiRankCoords _nRanksPerDir;
  }
  
  inline const MpiRankCoords &nRanksPerDir=resources::_nRanksPerDir;
  
}

#endif
