#ifndef _COMPONENTS_HPP
#define _COMPONENTS_HPP

#include "base/memory_manager.hpp"

namespace nissa
{
  /// Dynamic size
  constexpr Size DYNAMIC=-1;
  
  /// Specify the size at compile time
  template <Size SIZE=DYNAMIC>
  struct TensCompSize
  {
    /// Value beyond end
    static constexpr Size sizeAtCompileTime()
    {
      return SIZE;
    };
  };
  
  /// Row or column
  enum RwCl{RW,CL,ANY};
  
  /// Tensor component defined by base type S
  ///
  /// Inherit from S to get size
  template <typename S,
	    RwCl RC=RW,
	    int Which=0>
  struct TensCompIdx
  {
    /// Transposed type of component
    static constexpr RwCl TRANSP=(RC==ANY)?ANY:((RC==CL)?RW:CL);
    
    /// Transposed component
    using Transp=TensCompIdx<S,TRANSP,Which>;
    
    /// Base type
    typedef S Base;
    
    /// Value
    Size i;
    
    /// Check if the size is known at compile time
    static constexpr bool SizeIsKnownAtCompileTime=Base::sizeAtCompileTime()!=DYNAMIC;
    
    /// Init from int
    explicit constexpr TensCompIdx(Size i) : i(i)
    {
    }
    
    /// Default constructor
    TensCompIdx()
    {
    }
    
    /// Convert to actual value
    operator Size&()
    {
      return i;
    }
    
    /// Convert to actual value with const attribute
    operator const Size&() const
    {
      return i;
    }
    
    /// Transposed index
    auto transp() const
    {
      return Transp{i};
    }
  };
  
  /// Predicate returning whether the size is known ow not at compile time
  template <bool Asked=true>
  struct SizeIsKnownAtCompileTime
  {
    /// Internal implementation
    template <typename T>
    struct t
    {
      /// Predicate result
      static constexpr bool value=(T::SizeIsKnownAtCompileTime==Asked);
    };
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Define the base type for a component of name NAME and size SIZE
  #define DEFINE_BASE_COMP(NAME,SIZE)			\
  /*! Base type for a NAME index */			\
    struct _ ## NAME : public TensCompSize<SIZE>	\
  {							\
  }
  
  /// Provide the access to a component, through a function naed ACCESS
#define PROVIDE_COMP_ACCESS(ACCESS,...)		\
  auto inline ACCESS(const int64_t& i)		\
  {						\
    return __VA_ARGS__{i};			\
  }
  
  /// Define a component which cannot be included more than once
#define DEFINE_SINGLE_COMP_IDX(NAME,ACCESS,SIZE)	\
  DEFINE_BASE_COMP(NAME,SIZE);				\
  /*! NAME index */					\
  using NAME ## Idx=TensCompIdx<_ ## NAME,ANY,0>;	\
							\
  PROVIDE_COMP_ACCESS(ACCESS, NAME ## Idx)
  
  /// Define a component which can be included twice
#define DEFINE_RW_OR_CL_COMP_IDX(NAME,ACCESS,SIZE)	\
  DEFINE_BASE_COMP(NAME,SIZE);				\
  /*! NAME index */					\
  template <RwCl RC=RW,					\
	    int Which=0>				\
  using _ ## NAME ## Idx=TensCompIdx<_ ## NAME,RC,Which>;	\
							\
  using Rw ## NAME ## Idx = _ ## NAME ## Idx<RW,0>;	\
  using Cl ## NAME ## Idx = _ ## NAME ## Idx<CL,0>;	\
  using NAME ## Idx = Rw ## NAME ## Idx;	\
  							\
  PROVIDE_COMP_ACCESS(rw ## NAME, _ ## NAME ## Idx<RW,0>);	\
  PROVIDE_COMP_ACCESS(cl ## NAME, _ ## NAME ## Idx<CL,0>);	\
  PROVIDE_COMP_ACCESS(ACCESS,_ ## NAME ## Idx<RW,0>)
  
  /// Complex index
  DEFINE_SINGLE_COMP_IDX(Compl,complexAccess,2);
  DEFINE_RW_OR_CL_COMP_IDX(Color,cl,NCOL);
  DEFINE_RW_OR_CL_COMP_IDX(Spin,sp,NDIRAC);
  DEFINE_RW_OR_CL_COMP_IDX(Dir,dir,NDIM);
  DEFINE_SINGLE_COMP_IDX(LocVol,locVol,DYNAMIC);
  DEFINE_SINGLE_COMP_IDX(LocVolEvn,locVolEvn,DYNAMIC);
  DEFINE_SINGLE_COMP_IDX(LocVolOdd,locVolOdd,DYNAMIC);
  
  /// Collection of components
  template <typename...Tc>
  using TensComps=std::tuple<Tc...>;
  
  /// Real part access to complex
  constexpr ComplIdx re{0};
  
  /// Imaginary part access to complex
  constexpr ComplIdx im{1};
  
  /// Loop over a range
#define LOOP_RANGE(TYPE,NAME,MIN,MAX)					\
  TYPE NAME{MIN};NAME<MAX;NAME++
  
#define REIM(NAME) auto NAME : {re,im}
  
#define _ALL_RW_OR_CL_COLORS(NAME,RC) LOOP_RANGE(ColorIdx<RC>,NAME,0,NCL)
#define _ALL_RW_OR_CL_SPINS(NAME,RC) LOOP_RANGE(SpinIdx<RC>,NAME,0,NDIRAC)
#define _ALL_RW_OR_CL_DIRS(NAME,RC) LOOP_RANGE(DirIdx<RC>,NAME,0,NDIM)
  
#define ALL_RW_COLORS(NAME) _ALL_RW_OR_CL_COLORS(NAME,RW)
#define ALL_CL_COLORS(NAME) _ALL_RW_OR_CL_COLORS(NAME,CL)
#define ALL_COLORS(NAME) ALL_RW_COLORS(NAME)
  
#define ALL_RW_SPINS(NAME) _ALL_RW_OR_CL_SPINS(NAME,RW)
#define ALL_CL_SPINS(NAME) _ALL_RW_OR_CL_SPINS(NAME,CL)
#define ALL_SPINS(NAME) ALL_RW_SPINS(NAME)
  
#define ALL_RW_DIRS(NAME) _ALL_RW_OR_CL_DIRS(NAME,RW)
#define ALL_CL_DIRS(NAME) _ALL_RW_OR_CL_DIRS(NAME,CL)
#define ALL_DIRS(NAME) ALL_RW_DIRS(NAME)
  
#define ALL_LOC_SITES(NAME) LOOP_RANGE(LocVol,NAME,loc_vol)
#define ALL_EVN_LOC_SITES(NAME) LOOP_RANGE(LocVolEvn,NAME,loc_volh)
#define ALL_ODD_LOC_SITES(NAME) LOOP_RANGE(LocVolOdd,NAME,loc_volh)
}

#endif
