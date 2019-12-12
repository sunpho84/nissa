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
  
  /// Base type for a space index
  struct _Space : public TensCompSize<>
  {
  };
  
  /// Base type for a color index
  struct _Color : public TensCompSize<3>
  {
  };
  
  /// Base type for a spin index
  struct _Spin : public TensCompSize<4>
  {
  };
  
  /// Base type for a complex index
  struct _Compl : public TensCompSize<2>
  {
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
    explicit TensCompIdx(Size i) : i(i)
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
  
  /// Collection of components
  template <typename...Tc>
  using TensComps=std::tuple<Tc...>;
  
  /// Space index
  using SpaceIdx=TensCompIdx<_Space,ANY,0>;
  
  /// Spin index
  template <RowCol RC=ROW,
	    int Which=0>
  using SpinIdx=TensCompIdx<_Spin,RC,Which>;
  
  /// Color index
  template <RowCol RC=ROW,
	    int Which=0>
  using ColorIdx=TensCompIdx<_Color,RC,Which>;
  
  /// Complex index
  using ComplIdx=TensCompIdx<_Compl,ANY,0>;
}

#endif
