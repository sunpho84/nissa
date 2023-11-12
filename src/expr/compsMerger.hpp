#ifndef _COMPS_MERGER_HPP
#define _COMPS_MERGER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/compsMerger.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/exprRefOrVal.hpp>
#include <expr/mergedComps.hpp>
#include <expr/nodeDeclaration.hpp>
#include <expr/scalar.hpp>
#include <expr/subNodes.hpp>
#include <metaprogramming/detectableAs.hpp>
#include <metaprogramming/templateEnabler.hpp>
#include <routines/ios.hpp>
#include <tuples/tupleFilter.hpp>
#include <tuples/tupleHasType.hpp>

namespace nissa
{
  PROVIDE_DETECTABLE_AS(CompsMerger);
  
  /// Components merger
  ///
  /// Forward declaration to capture the components
  template <typename _MC,
	    typename _E,
	    typename _Comps,
	    typename _Fund>
  struct CompsMerger;
  
#define THIS					\
  CompsMerger<CompsList<Mc...>,std::tuple<_E...>,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Components merger
  ///
  template <typename...Mc,
	    typename..._E,
	    typename...C,
	    typename _Fund>
  struct THIS :
    DetectableAsCompsMerger,
    SubNodes<_E...>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    static_assert(sizeof...(_E)==1,"Expecting 1 argument");
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    using ThisMergedComp=
      MergedComp<CompsList<Mc...>>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    IMPORT_SUBNODE_TYPES;
    
    /// Merged type
    using MergedExpr=SubNode<0>;
    
#define PROVIDE_MERGED_EXPR(ATTRIB)			\
    /*! Returns the merged expression */		\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
    ATTRIB MergedExpr& mergedExpr() ATTRIB		\
    {							\
      return SUBNODE(0);				\
    }
    
    PROVIDE_MERGED_EXPR(/* non const */);
    
    PROVIDE_MERGED_EXPR(const);
    
#undef PROVIDE_MERGED_EXPR
    
    static constexpr bool mergedCompHasDynamicSize=
      ThisMergedComp::sizeIsKnownAtCompileTime;
    
    /// Extra dynamic comp w.r.t merged expression
    using ExtraDynamicComps=
      std::conditional_t<mergedCompHasDynamicSize,
      CompsList<ThisMergedComp>,
      CompsList<>>;
    
    ExtraDynamicComps extraDynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return tupleGetSubset<typename CompsMerger::DynamicComps>(std::tuple_cat(mergedExpr().getDynamicSizes(),extraDynamicSizes));
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    bool canAssign()
    {
      return mergedExpr().canAssign();
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    This& operator=(const This& oth)
    {
      *static_cast<Base*>(this)=*static_cast<const Base*>(&oth);
      
      return *this;
    }
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      MergedExpr::canAssignAtCompileTime;
    
//     /////////////////////////////////////////////////////////////////
    
#define PROVIDE_RECREATE_FROM_EXPR(ATTRIB)			\
    /*! Returns a ATTRIB similar version */			\
    template <typename T>					\
    INLINE_FUNCTION						\
    decltype(auto) recreateFromExprs(T&& t) ATTRIB		\
    {								\
      return mergeComps<Mc...>(std::forward<T>(t));		\
    }
    
    PROVIDE_RECREATE_FROM_EXPR(/* non const */);
    
    PROVIDE_RECREATE_FROM_EXPR(const);
    
#undef PROVIDE_RECREATE_FROM_EXPR
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return recreateFromExprs(mergedExpr().getRef());		\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_EVAL(ATTRIB)						\
    template <typename...U>						\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION			\
    decltype(auto) eval(const U&...cs) ATTRIB				\
    {									\
      const auto procComp=						\
	[this]<typename Ui>(const Ui& c)				\
	{								\
	  if constexpr(std::is_same_v<Ui,ThisMergedComp>)		\
	    return c.decompose(mergedExpr().getDynamicSizes());		\
	  else								\
	    return std::make_tuple(c);					\
	};								\
      									\
      return								\
      std::apply(mergedExpr().eval,std::tuple_cat(procComp(cs)...));	\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/*non const*/);
    
#undef PROVIDE_EVAL
    
    /// Construct
    template <typename T>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    CompsMerger(T&& arg) :
      SubNodes<_E...>(std::forward<T>(arg))
    {
      if constexpr(mergedCompHasDynamicSize)
	std::get<ThisMergedComp>(extraDynamicSizes)=mergedExpr().template getMergedCompsSize<ThisMergedComp>();
    }
  };
  
  /// Merges a subset of components
  template <typename MC,
	    typename _E,
	    ENABLE_THIS_TEMPLATE_IF(isNode<_E>)>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  auto mergeComps(_E&& e)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    /// Type returned when evaluating the expression
    using Fund=
      typename E::Fund;
    
    /// Visible components
    using Comps=
      CompsMerge<MC,typename E::Comps>;
    
    return
      CompsMerger<MC,
		  std::tuple<decltype(e)>,
		  Comps,
		  Fund>(std::forward<_E>(e));
  }
}

#endif
