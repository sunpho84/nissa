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
    template <DerivedFromComp...D,
	      DerivedFromComp...C>
    static INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    MergedComp merge(const std::tuple<D...>& dynamicSizes,
		     const C&...c)
    {
      return orderedIndex<Cp...>(dynamicSizes,c...);
    }
    
    /// Gets the components of a merged component
    template <typename D>
    static INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    Comps decompose(const D& dynamicSizes,
		    const MergedComp& i)
    {
      return indexDecompose<Comps>(dynamicSizes,i());
    }
    
    /// Gets the components of this merged component
    template <typename D>
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    Comps decompose(const D& dynamicSizes) const
    {
      return MergedComp::decompose(dynamicSizes,*this);
    }
    
    /// Gets the components of this merged component
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    Comps decompose() const
    {
      return this->decompose(std::make_tuple(),*this);
    }
    
    /// Initialize from dynamicSizes and components
    template <DerivedFromComp...D,
	      DerivedFromComp...E>
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    MergedComp(const std::tuple<D...>& dynamicSizes,
	       const E&...e) :
      MergedComp(MergedComp::merge(dynamicSizes,e...))
    {
    }
    
    /// Initialize from components
    template <DerivedFromComp...E>
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    MergedComp(const E&...e) :
      MergedComp(std::tuple<>{},e...)
    {
    }
    
    /// Default copy constructor
    INLINE_FUNCTION constexpr
    MergedComp(const MergedComp&)=default;
  };
  
  template <typename Tp,
	    typename T>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  MergedComp<Tp> mergedComp(T&& i)
  {
    return i;
  };
  
  /////////////////////////////////////////////////////////////////
  
  namespace impl
  {
    /// Merge components
    ///
    /// Forward declaration
    template <typename PastComps,
	      typename ToMerge,
	      typename CompsToProcess,
	      bool MergedHasAlreadyBeenInserted>
    struct _CompsMerge;
    
    /// Merge components
    ///
    /// Case in which the processed component matches the first of the
    /// components to be merged
    template <typename...PastComps,
	      typename...ToMerge,
	      typename HeadToProcess,
	      typename...TailToProcess,
	      bool MergedHasAlreadyBeenInserted>
    struct _CompsMerge<CompsList<PastComps...>,
		       CompsList<ToMerge...>,
		       CompsList<HeadToProcess,TailToProcess...>,
		       MergedHasAlreadyBeenInserted>
    {
      static constexpr bool hasHeadToProcessToBeMerged=
	(std::is_same_v<HeadToProcess,ToMerge> or...);
      
      /// Resulting components: if the merged components have been
      /// inserted, drops the first component to be merged, otherwise
      /// insert all the components
      using ResComps=
	std::conditional_t<hasHeadToProcessToBeMerged,
			   std::conditional_t<MergedHasAlreadyBeenInserted,
					      CompsList<PastComps...>,
					      CompsList<PastComps...,MergedComp<CompsList<ToMerge...>>>>,
			   CompsList<PastComps...,HeadToProcess>>;
      
      /// Resulting type, detected iteratively, pointing out that the
      /// merged component has for sure been inserted
      using type=
	typename _CompsMerge<ResComps,
			     CompsList<ToMerge...>,
			     CompsList<TailToProcess...>,
			     MergedHasAlreadyBeenInserted or hasHeadToProcessToBeMerged>::type;
    };
    
    /// Merge components
    ///
    /// Case in which no component needs to be processed
    template <typename C,
	      typename M>
    struct _CompsMerge<C,
		       M,
		       CompsList<>,
		       true>
    {
      /// Resulting type obtained passing all components to be
      /// processed to the output
      using type=C;
    };
  }
  
  /// Merge compontents
  template <typename MC,
	    typename C>
  using CompsMerge=
    impl::_CompsMerge<CompsList<>,MC,C,false>::type;
}

#endif
