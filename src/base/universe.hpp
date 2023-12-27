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
  
  DECLARE_TRANSPOSABLE_COMP(PerpDir,int,NDIM-1,perpDir);
  
  DECLARE_UNTRANSPOSABLE_COMP(Ori,int,2,createOri);
  
  /// Coordinates for a given type
  template <typename F>
  using Coords=
    StackTens<CompsList<Dir>,F>;
  
  /// Index of time direction
  static constexpr Dir timeDir=0;
  
  /// Mapping of scidac to nissa directions
  CUDA_DEVICE static constexpr Coords<Dir> scidacNissaDirMapping=
    [](const Dir& in)
    {
      if(in==timeDir)
	return in;
      else
	return NDIM-in;
    };
  
  /// Decompose lx index f into the coordinates in sizes s
  template <typename F,
	    typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
  Coords<C> decomposeLxToCoords(F f,
				const Coords<C>& _s)
  {
    Coords<C> res;
    const Coords<F>& s=_s.template reinterpretFund<F>();
    Coords<F>& work=
      res.template reinterpretFund<F>();
    
    for(Dir mu=NDIM-1;mu>=0;mu--)
      {
	F t=s[mu];
	work[mu]=f%t;
	f/=t;
      }
    
    return res;
  }
  
  /// Finds the index f of coordinates in sizes s
  template <typename F,
	    typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
  F lxOfCoords(const Coords<C>& _x,
	       const Coords<C>& _s)
  {
    const Coords<F>& x=_x.template reinterpretFund<F>();
    const Coords<F>& s=_s.template reinterpretFund<F>();
    
    F f=0;
    
    for(Dir mu=0;mu<NDIM;mu++)
      f=f*s[mu]+x[mu];
    
    return f;
  }
  
  /// Backward, see real imag comment
#define bw Ori(0)
  
  /// Forward
#define fw Ori(1)
  
  /// Number of dimensions
#define nDim Dir(NDIM)
  
  /// Perpendicular direction of a given dir
  constexpr StackTens<CompsList<Dir,PerpDir>,Dir> perpDirOf=
    [] (const Dir& dir,
	const PerpDir& perpDir) -> Dir
    {
      return (perpDir()>=dir())?(perpDir()+1):perpDir();
    };
}

#endif
