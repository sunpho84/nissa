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
#include <expr/subExprs.hpp>
#include <routines/ios.hpp>
#include <tuples/tupleFilter.hpp>
#include <tuples/tupleHasType.hpp>

namespace nissa
{
  PROVIDE_FEATURE(CompsMerger);
  
  /// Components merger
  ///
  /// Forward declaration to capture the components
  template <typename _MC,
	    typename _E,
	    typename _Comps,
	    typename _Fund>
  struct CompsMerger;
  
#define THIS					\
  CompsMerger<CompsList<Mc...>,_E,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Components merger
  ///
  template <typename...Mc,
	    typename _E,
	    typename...C,
	    typename _Fund>
  struct THIS :
    CompsMergerFeat,
    SingleSubExpr<THIS>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    using ThisMergedComp=
      MergedComp<CompsList<Mc...>>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    /// Expression to merge
    NodeRefOrVal<_E> subExpr;
    
    /// Type of the merged expr
    using MergedExpr=std::decay_t<_E>;
    
    /// Exec as the merged expression
    static constexpr ExecSpace execSpace=
      MergedExpr::execSpace;
    
    static constexpr bool mergedCompHasDynamicSize=
      not ThisMergedComp::sizeIsKnownAtCompileTime;
    
    /// Extra dynamic comp w.r.t merged expression
    using ExtraDynamicComps=
      std::conditional_t<mergedCompHasDynamicSize,
      CompsList<ThisMergedComp>,
      CompsList<>>;
    
    ExtraDynamicComps extraDynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    decltype(auto) getDynamicSizes() const
    {
      return tupleGetSubset<typename CompsMerger::DynamicComps>(std::tuple_cat(subExpr.getDynamicSizes(),extraDynamicSizes));
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    bool canAssign()
    {
      return subExpr.canAssign();
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
    
    /////////////////////////////////////////////////////////////////
    
    /// Returns a similar version
    template <typename T>
    INLINE_FUNCTION constexpr
    auto recreateFromExprs(T&& t) const
    {
      return mergeComps<Mc...>(std::forward<T>(t));
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return recreateFromExprs(subExpr.getRef());		\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      CompsMerger<CompsList<Mc...>,
		  typename _E::template ReinterpretFund<NFund>,
		  CompsList<C...>,
		  NFund>;
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_EVAL(ATTRIB)						\
    template <typename...U>						\
    HOST_DEVICE_ATTRIB constexpr INLINE_FUNCTION			\
    decltype(auto) eval(const U&...cs) ATTRIB				\
    {									\
      const auto procComp=						\
	[this]<typename Ui>(const Ui& c)				\
	{								\
	  if constexpr(std::is_same_v<Ui,ThisMergedComp>)		\
	    return c.decompose(subExpr.getDynamicSizes());		\
	  else								\
	    return std::make_tuple(c);					\
	};								\
									\
      return								\
	std::apply([this](const auto&...c) ->decltype(auto)		\
	{								\
	  return subExpr.eval(c...);				\
	},std::tuple_cat(procComp(cs)...));				\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/*non const*/);
    
#undef PROVIDE_EVAL
    
    /// Construct
    template <typename T>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    CompsMerger(T&& arg)
      requires(std::is_same_v<std::decay_t<T>,std::decay_t<_E>>)
      : subExpr{std::forward<T>(arg)}
    {
      if constexpr(mergedCompHasDynamicSize)
	std::get<ThisMergedComp>(extraDynamicSizes)=subExpr.template getMergedCompsSize<CompsList<Mc...>>();
    }
  };
  
  /// Merges a subset of components
  template <typename MC,
	    DerivedFromNode _E>
  HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
  decltype(auto) mergeComps(_E&& e)
  {
    if constexpr(std::tuple_size_v<MC> <=1)
      return (NodeRefOrVal<_E>)e;
    else
      {
	/// Base passed type
	using E=
	  std::decay_t<_E>;
	
#ifndef COMPILING_FOR_DEVICE
	master_printf("getting a lvalue: %d rvalue: %d\n",std::is_lvalue_reference_v<decltype(e)>,std::is_rvalue_reference_v<decltype(e)>);
#endif
	
	/// Type returned when evaluating the expression
	using Fund=
	  typename E::Fund;
	
	/// Visible components
	using Comps=
	  CompsMerge<MC,typename E::Comps>;
	
	return
	  CompsMerger<MC,
		      decltype(e),
		      Comps,
		      Fund>(std::forward<_E>(e));
      }
  }
}

#endif
