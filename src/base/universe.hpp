#ifndef _UNIVERSE_HPP
#define _UNIVERSE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <expr/comps.hpp>
#include <expr/stackTens.hpp>

namespace nissa
{
  DECLARE_UNTRANSPOSABLE_COMP(Parity,int,2,createParity);
  DECLARE_TRANSPOSABLE_COMP(Dir,int,NDIM,dir);
  
  DECLARE_UNTRANSPOSABLE_COMP(Ori,int,2,createOri);
  
  /// Coordinates for a given type
  template <DerivedFromComp Comp>
  using Coords=
    StackTens<CompsList<Dir>,Comp>;
  
  /// Backward, see real imag comment
#define bw Ori(0)
  
  /// Forward
#define fw Ori(1)
  
  /// Number of dimensions
#define nDim Dir(NDIM)
}

#endif
