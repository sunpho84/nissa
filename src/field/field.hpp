#ifndef _FIELD_HPP
#define _FIELD_HPP

#include "component.hpp"
#include "tensor.hpp"

namespace nissa
{
  ///////////// TO BE MOVED SOMEWHERE ELSE ////////
  template <TensStorageLocation>
  struct _ParallelLoop;
  
  template <>
  struct _ParallelLoop<TensStorageLocation::ON_CPU>
  {
    template <typename S,
	      typename F>
    static void loop(const S& min,const S& max,F&& f)
    {
      for(S i=min;i<max;i++)
	f(i);
    }
  };
  
#ifdef USE_CUDA
  template <>
  struct _ParallelLoop<TensStorageLocation::ON_GPU>
  {
    template <typename S,
	      typename F>
    static void loop(const S& min,const S& max,F&& f)
    {
      cuda_parallel_for(__LINE__,__FILE__,min,max,f);
    }
  };
#endif
  
  template <TensStorageLocation SL,
	    typename S,
	    typename F>
  void parallelLoop(const S& min,const S& max,F&& f)
  {
    _ParallelLoop<SL>::loop(min,max,std::forward<F>(f));
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Ordinary storage location used when not specified otherwise
  static constexpr TensStorageLocation OrdinaryFieldStorageLocation=
#ifdef USE_CUDA
    TensStorageLocation::ON_GPU
#else
    TensStorageLocation::ON_CPU
#endif
    ;
  
  /// Internal layout specifier
  enum class FieldLayout{CPU,GPU,VECT};
  
  /// Specify which kind of halo to use
  enum class HaloKind{NO_HALO,WITH_BORDERS,WITH_EDGES};
  
  /// Properties of the SpaceTime component
  template <typename T>
  struct SpaceTimeProps
  {
    inline static T locSize();
    inline static T bordSize();
    inline static T edgeSize();
    
    inline static T totSize(HaloKind haloKind)
    {
      Size res=0;
      
      switch(haloKind)
	{
	case HaloKind::NO_HALO:
	  res=locSize();
	  break;
	case HaloKind::WITH_BORDERS:
	  res=locSize()+bordSize();
	  break;
	case HaloKind::WITH_EDGES:
	  res=locSize()+bordSize()+edgeSize();
	  break;
	}
      
      return static_cast<T>(res);
    }
  };
  
  /// Define the local properties of the given SpaceTime kind
#define PROVIDE_LOC_PROP(TYPE,LOC_SIZE,BORD_SIZE,EDGE_SIZE)	\
  template <>							\
  inline TYPE SpaceTimeProps<TYPE>::locSize()			\
  {								\
    return TYPE{LOC_SIZE};					\
  }								\
								\
  template <>							\
  inline TYPE SpaceTimeProps<TYPE>::bordSize()			\
  {								\
    return TYPE{BORD_SIZE};					\
  }								\
								\
								\
  template <>							\
  inline TYPE SpaceTimeProps<TYPE>::edgeSize()			\
  {								\
    return TYPE{EDGE_SIZE};					\
  }
  
  PROVIDE_LOC_PROP(LocVolIdx,loc_vol,bord_vol,edge_vol)
  PROVIDE_LOC_PROP(LocVolEvnIdx,loc_volh,bord_volh,edge_volh)
  PROVIDE_LOC_PROP(LocVolOddIdx,loc_volh,bord_volh,edge_volh)
  
#undef PROVIDE_LOC_PROP
  
  /// Ordinary layout used when not specified otherwise
  static constexpr FieldLayout OrdinaryLayout=
#ifdef USE_CUDA
    FieldLayout::GPU
#else
    FieldLayout::CPU
#endif
    ;
  
  /// Contains the field halo properties
  template <HaloKind>
  struct HaloFieldProperties;
  
  /// Field type
  ///
  /// Internal layout is delegated to CompsTraits
  template <typename _SPT,                      // Spacetime component
	    typename _Tc,                       // Other components
	    typename _F=double,                 // Fundamental type
	    FieldLayout FL=OrdinaryLayout,      // Layout
	    TensStorageLocation SL=OrdinaryFieldStorageLocation,
	    bool IsRef=false>
  struct Field : public Subscribable<Field<_SPT,_Tc,_F,FL,SL,IsRef>>
  {
    /// Components traits
    using CT=CompsTraits<Field<_SPT,_Tc,_F,FL,SL,IsRef>>;
    
    /// Kind of halo
    const HaloKind haloKind;
    
    /// Storing data
    Tens<typename CT::Comps,typename CT::F,SL,IsRef> data;
    
    template <typename F>
    static void spaceTimeParallelLoop(F&& f)
    {
      using S=typename CT::SpaceTimeComp;
      
      parallelLoop<SL>((S)0,SpaceTimeProps<S>::locSize(),std::forward<F>(f));
    }
    
    /// Create taking the dynamical sizes as argument
    template <typename...D>
    Field(D&&...d) : Field(HaloKind::NO_HALO,std::forward<D>(d)...)
    {
    }
    
    /// Create taking the dynamical sizes as argument
    template <typename...D>
    Field(HaloKind haloKind,
	  D&&...d) :
      haloKind(haloKind),
      data(SpaceTimeProps<typename CT::SpaceTimeComp>::totSize(haloKind),std::forward<D>(d)...)
    {
    }
    
    /// Evaluate
    template <typename...A>
    const auto& eval(A&&...a) const
    {
      return data(std::forward<A>(a)...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(eval);
  };
  
  /// Field inner types
  ///
  /// Forward definition
  template <FieldLayout FL,
	    typename _S,
	    typename _F,
	    typename _Tc>
  struct FieldInnerTypes;
  
  /// Field inner types
  ///
  /// CPU case
  template <typename _S,
	    typename _F,
	    typename...Tp>
  struct FieldInnerTypes<FieldLayout::CPU,_S,_F,TensComps<Tp...>>
  {
    /// Spacetime: no change
    using S=_S;
    
    /// Fundamental type: no change
    using F=_F;
    
    /// Components: put spacetime before all the rest
    using Comps=TensComps<S,Tp...>;
  };
  
  /// Field inner types
  ///
  /// GPU case
  template <typename _S,
	    typename _F,
	    typename...Tp>
  struct FieldInnerTypes<FieldLayout::GPU,_S,_F,TensComps<Tp...>>
  {
    /// Spacetime: no change
    using S=_S;
    
    /// Fundamental type: no change
    using F=_F;
    
    /// Components: put spacetime after all the rest
    using Comps=TensComps<Tp...,S>;
  };
  
  /// Properties of components of field
  template <typename _S,
	    typename _Tc,
	    typename _F,
	    FieldLayout FL,
	    TensStorageLocation SL>
  struct CompsTraits<Field<_S,_Tc,_F,FL,SL>>
  {
    /// Actual type
    using Type=Field<_S,_Tc,_F,FL,SL>;
    
    /// List of valid types
    using SpaceTimeTypes=std::tuple<LocVolIdx,LocVolEvnIdx,LocVolOddIdx>;
    
    static_assert(std::tuple_size<TupleCommonTypes<std::tuple<_S>,SpaceTimeTypes>>::value==1,"Invalid spacetime component");
    
    /// Field inner types
    using FIT=FieldInnerTypes<FL,_S,_F,_Tc>;
    
    /// Spacetime inner type
    using SpaceTimeComp=typename FIT::S;
    
    /// Components
    using Comps=typename FIT::Comps;
    
    /// Fundamental type
    using F=typename FIT::F;
  };
}

#endif
