#ifndef _INDEXCOMPUTER_HPP
#define _INDEXCOMPUTER_HPP

#include <tensor/component.hpp>

namespace nissa
{
  template <typename Comps>
  struct IndexComputer;
  
  template <typename...TC>
  struct IndexComputer<TensorComps<TC...>>
  {
    /// Type to be used for the index
    using Index=
      std::common_type_t<int,typename TC::Index...>;
    
    /// List of all statically allocated components
    using StaticComps=
      TupleFilter<predicate::SizeIsKnownAtCompileTime<true>::t,TensorComps<TC...>>;
    
    /// List of all dynamically allocated components
    using DynamicComps=
      TupleFilter<predicate::SizeIsKnownAtCompileTime<false>::t,UniqueTuple<TC...>>;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes{};
    
    /// Set the dynamic sizes
    template <typename...TD>
    CUDA_HOST_DEVICE
    void setDynamicSizes(const TensorComps<TD...>& td)
    {
      static_assert(sizeof...(TD)==std::tuple_size_v<DynamicComps>,"Cannot allocate without knowledge of all the dynamic sizes");
      
      tupleFillWithSubset(dynamicSizes,td);
    }
    
    /// Static component maximal value
    static constexpr Index staticPartMaxValue=
      ((TC::sizeIsKnownAtCompileTime?
	TC::sizeAtCompileTime():
       Index{1})*...);
    
    /// Size of the Tv component
    ///
    /// Case in which the component size is known at compile time
    template <typename Tv,
	      ENABLE_THIS_TEMPLATE_IF(Tv::sizeIsKnownAtCompileTime)>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    constexpr const typename Tv::Index compSize()
      const
    {
      return
	Tv::sizeAtCompileTime();
    }
    
    /// Size of the Tv component
    ///
    /// Case in which the component size is not knwon at compile time
    template <typename Tv,
	      ENABLE_THIS_TEMPLATE_IF(not Tv::sizeIsKnownAtCompileTime)>
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    const typename Tv::Index& compSize()
      const
    {
      return std::get<Tv>(dynamicSizes)();
    }
    
    /// Calculate the index - no more components to parse
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    const Index& _index(const Index& outer) ///< Value of all the outer components
      const
    {
      return outer;
    }
    
    /// Determine whether the components are all static, or not
    static constexpr bool allCompsAreStatic=
      std::is_same<DynamicComps,std::tuple<>>::value;
    
    // /// Computes the maximal value size at compile time, if known
    // static constexpr Index maxValAtCompileTime=
    //   allCompsAreStatic?(Index)staticPartMaxValue:(Index)DYNAMIC;
    
    /// Compute the maximal value at compile time
    constexpr Index maxVal()
    {
      return ((compSize<TC>())*...);
    }
    
    /// Calculate index iteratively
    ///
    /// Given the components (i,j,k) we must compute ((0*ni+i)*nj+j)*nk+k
    ///
    /// The parsing of the variadic components is done left to right, so
    /// to compute the nested bracket list we must proceed inward. Let's
    /// say we are at component j. We define outer=(0*ni+i) the result
    /// of inner step. We need to compute thisVal=outer*nj+j and pass it
    /// to the inner step, which incorporate iteratively all the inner
    /// components. The first step requires outer=0.
    template <typename T,
    	      typename...Tp>
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    Index _index(const Index& outer, ///< Value of all the outer components
			    T&& thisComp,       ///< Currently parsed component
			    Tp&&...innerComps)  ///< Inner components
      const
    {
      /// Remove reference and all attributes to access to types
      using Tv=
	std::decay_t<T>;
      
      /// Size of this component
      const Index thisSize=
	compSize<Tv>();
      
      /// Value of the index when including this component
      const Index thisVal=
	outer*thisSize+thisComp();
      
      return
	_index(thisVal,innerComps...);
    }
    
    /// Dispatch the internal index calculation
    ///
    /// This works when the passed components are already well ordered
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    Index operator()(const TC&...comps) const
    {
      return _index(Index(0),comps...);
    }
  };
}

#endif
