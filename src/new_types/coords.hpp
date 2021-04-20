#ifndef _COORDS_HPP
#define _COORDS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <tensor/tensor.hpp>

namespace nissa
{
  DECLARE_COMPONENT(Dir,int,NDIM);
  
  /// Temporal direction
  inline CUDA_DEVICE Dir tDir=0;
  
  /// X direction
  inline CUDA_DEVICE Dir xDir=1;
  
  /// Y direction
  inline CUDA_DEVICE Dir yDir= 2;
  
  /// Z direction
  inline CUDA_DEVICE Dir zDir=3;
  
#define FOR_ALL_DIRS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(Dir,NAME)
  
#define FOR_ALL_SPATIAL_DIRS(NAME)	\
  for(Dir NAME=1;NAME<NDIM;NAME++)
  
  template <typename T>
  using Coords=Tensor<OfComps<Dir>,T>;
  
  using Momentum=Tensor<OfComps<Dir>,double>;
  
  /// Copy coordinates
  template <typename T>
  void coord_copy(Coords<T>& out,const Coords<T>& in)
  {
    FOR_ALL_DIRS(mu)
      out(mu)=in(mu);
  }
  
  /// Summ the coordinates a1 and a2 modulo l
  template <typename T>
  void coord_summ(Coords<T>& s,const Coords<T>& a1,const Coords<T>& a2,const Coords<T>& l)
  {
    FOR_ALL_DIRS(mu)
      s(mu)=(a1(mu)+a2(mu))%l(mu);
  }
  
  /// Summassgn the coordinates s and a, modulo l
  template <typename T>
  void coord_summassign(Coords<T>& s,const Coords<T>& a,const Coords<T>& l)
  {
    coord_summ(s,s,a,l);
  }
}

#endif
