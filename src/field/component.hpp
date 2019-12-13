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
  enum RowCol{ROW,COL,ANY};
  
  /// Tensor component defined by base type S
  ///
  /// Inherit from S to get size
  template <typename S,
	    RowCol RC=ROW,
	    int Which=0>
  struct TensCompIdx
  {
    /// Transposed type of component
    static constexpr RowCol TRANSP=(RC==ANY)?ANY:((RC==COL)?ROW:COL);
    
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
    
    /// Convert to actual value
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
  
  /// Provide the access to a component
#define PROVIDE_COMP_ACCESS(ACCESS,...)		\
  auto ACCESS(const int64_t& i)			\
  {						\
    return __VA_ARGS__{i};			\
  }
  
  /// Define a component which cannot be included more than once
#define DEFINE_SINGLE_COMP_IDX(NAME,ACCESS,SIZE)	\
  DEFINE_BASE_COMP(NAME,SIZE);				\
  /*! NAME index */					\
  using NAME ## Idx=TensCompIdx<_ ## NAME,ANY,0>;	\
							\
  PROVIDE_COMP_ACCESS(ACCESS,NAME ## Idx)
  
#define DEFINE_ROW_OR_COL_COMP_IDX(NAME,ACCESS,SIZE)	\
  DEFINE_BASE_COMP(NAME,SIZE);				\
  /*! NAME index */					\
  template <RowCol RC=ROW,				\
	    int Which=0>				\
  using NAME ## Idx=TensCompIdx<_ ## NAME,RC,Which>;	\
							\
  PROVIDE_COMP_ACCESS(row ## NAME,NAME ## Idx<ROW,0>);	\
  PROVIDE_COMP_ACCESS(col ## NAME,NAME ## Idx<COL,0>);	\
  PROVIDE_COMP_ACCESS(ACCESS,NAME ## Idx<ROW,0>)
  
  /// Complex index
  DEFINE_SINGLE_COMP_IDX(Compl,complexAccess,2);
  DEFINE_ROW_OR_COL_COMP_IDX(Color,cl,NCOL);
  DEFINE_ROW_OR_COL_COMP_IDX(Spin,sp,NDIRAC);
  DEFINE_ROW_OR_COL_COMP_IDX(Lorentz,lorentz,NDIM);
  DEFINE_SINGLE_COMP_IDX(LocVol,locVol,DYNAMIC);
  DEFINE_SINGLE_COMP_IDX(LocVolEvn,locVolEvn,DYNAMIC);
  DEFINE_SINGLE_COMP_IDX(LocVolOdd,locVolOdd,DYNAMIC);
  
  /// Collection of components
  template <typename...Tc>
  using TensComps=std::tuple<Tc...>;
}

#endif
