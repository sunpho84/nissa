#ifndef _COORDS_HPP
#define _COORDS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <tensor/tensor.hpp>

namespace nissa
{
  DECLARE_COMPONENT(Orientation,int,2);
  
  /// Bakcward direction
#define BACW Orientation(0)
  
  /// Forward direction
#define FORW Orientation(1)
    
  DECLARE_COMPONENT(Dir,int,NDIM);
  
  DECLARE_COMPONENT(PerpDir,int,NDIM-1);
  
  /// Convenient alias for Dir
  using Lorentz=
    Dir;
  
  /// Temporal direction
#define tDir Dir(0)
  
  /// X direction
#define xDir Dir(1)
  
  /// Y direction
#define yDir Dir(2)
  
  /// Z direction
#define zDir Dir(3)
  
#define FOR_ALL_DIRS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(Dir,NAME)
  
#define FOR_ALL_PERP_DIRS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(PerpDir,NAME)
  
#define UNROLL_FOR_ALL_DIRS(NAME,CORE...)	\
  FOR_ALL_COMPONENT_VALUES(Dir,NAME,CORE)
  
#define FOR_ALL_SPATIAL_DIRS(NAME)				\
  FOR_ALL_COMPONENT_VALUES_STARTING_AT(Dir,NAME,1)
  
#define UNROLL_FOR_ALL_SPATIAL_DIRS(NAME,CORE...)				\
  UNROLL_FOR_ALL_COMPONENT_VALUES_STARTING_AT(Dir,NAME,1,CORE)
  
   /// Coordinates
  template <typename T,
	    StorLoc SL=DefaultStorage>
  using Coords=Tensor<OfComps<Dir>,T,SL,TensorDynamicity::STACKED_TENSOR>;
  
  /// Momementum
  using Momentum=Tensor<OfComps<Dir>,double,DefaultStorage,TensorDynamicity::STACKED_TENSOR>;
  
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
  
  /////////////////////////////////////////////////////////////////
  
  /// Clifford algebra elements nasty remove from here
#define NGAMMA 15
  
  DECLARE_COMPONENT(Gamma,int,NGAMMA);
  
  inline CUDA_DEVICE Gamma GammaId=0;
  
  inline CUDA_DEVICE Gamma GammaX=1;
  
  inline CUDA_DEVICE Gamma GammaY=2;
  
  inline CUDA_DEVICE Gamma GammaZ=3;
  
  inline CUDA_DEVICE Gamma GammaT=4;
  
  inline CUDA_DEVICE Gamma Gamma5=5;
}

#endif
