#ifndef _MERGEDCOMPS_HPP
#define _MERGEDCOMPS_HPP

#include <expr/comps.hpp>
#include <expr/indexComputer.hpp>

namespace nissa
{
  /////////////////////////////////////////////////////////////////
  
  /// Component obtained from merging other components
  ///
  /// Forward declaration
  template <typename T>
  struct MergedComp;
  
#define THIS MergedComp<CompsList<Cp...>>
  
#define BASE BaseComp<THIS,decltype(((Cp{}())*...*1)),(Cp::sizeAtCompileTime*...*1)>
  
  /// Component obtained from merging other components
  template <typename...Cp>
  struct THIS :
    BASE
  {
    static_assert((isComp<Cp> and ...),"Cannot merge other than components");
    
    using Base=BASE;
    
#undef BASE
    
#undef THIS
    
    using Base::Base;
    
    using Comps=
      CompsList<Cp...>;
    
    /// Returns the merged component from the unmerged one
    template <typename D,
	      typename...E>
    static INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    MergedComp merge(const D& dynamicSizes,
		     const CompFeat<E>&...e)
    {
      return orderedIndex<Comps>(dynamicSizes,~e...);
    }
    
    /// Gets the components of a merged component
    template <typename D>
    static INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Comps decompose(const D& dynamicSizes,
		    const MergedComp& i)
    {
      return indexDecompose<Comps>(dynamicSizes,i());
    }
    
    /// Gets the components of this merged component
    template <typename D>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Comps decompose(const D& dynamicSizes) const
    {
      return MergedComp::decompose(dynamicSizes,*this);
    }
    
    /// Gets the components of this merged component
    template <typename D>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Comps decompose() const
    {
      return this->decompose(std::make_tuple(),*this);
    }
  };
  
  template <typename Tp,
	    typename T>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  MergedComp<Tp> mergedComp(T&& i)
  {
    return i;
  };
}

#endif
