#ifndef _COORDS_HPP
#define _COORDS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <tensor/tensor.hpp>

namespace nissa
{
  DECLARE_COMPONENT(Direction,int,NDIM);
  
  /// Temporal direction
  constexpr Direction timeDirection=0;
  
  //DECLARE_TEMPLATED_COMPONENT(Coords,NDIM);
  
#define FOR_ALL_DIRECTIONS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(Direction,NAME)
  
  template <typename T>
  using Coords=Tensor<OfComps<Direction>,T>;
  
  using Momentum=Tensor<OfComps<Direction>,double>;
  
  /// Copy coordinates
  template <typename T>
  void coord_copy(Coords<T>& out,const Coords<T>& in)
  {
    FOR_ALL_DIRECTIONS(mu)
      out(mu)=in(mu);
  }
  
  /// Summ the coordinates a1 and a2 modulo l
  template <typename T>
  void coord_summ(Coords<T>& s,const Coords<T>& a1,const Coords<T>& a2,const Coords<T>& l)
  {
    FOR_ALL_DIRECTIONS(mu)
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
