#ifndef _COMPREDUCE_HPP
#define _COMPREDUCE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/compReducer.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/node.hpp>
#include <metaprogramming/arithmeticOperatorsViaCast.hpp>
#include <tuples/tupleCat.hpp>

namespace nissa
{
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(CompReducer);
  
  /// CompReducer
  ///
  /// Forward declaration to capture the components
  template <DerivedFromComp Rc,
	    DerivedFromNode _E,
	    typename _Comps,
	    typename _Fund,
	    typename Combiner>
  struct CompReducer;
  
#define THIS					\
  CompReducer<Rc,_E,CompsList<C...>,_Fund,Combiner>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// CompReducer
  ///
  template <DerivedFromComp Rc,
	    DerivedFromNode _E,
	    DerivedFromComp...C,
	    typename _Fund,
	    typename Combiner>
  struct THIS :
    CompReducerFeat,
    SingleSubExpr<THIS>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Components
    using ReduceComp=
      Rc;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// Expression to reduce
    NodeRefOrVal<_E> subExpr;
    
    /// Type of the component reduced expression
    using CompReducedExpr=std::decay_t<_E>;
    
    /// Executes where does the reference
    static constexpr ExecSpace execSpace=
      CompReducedExpr::execSpace;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    decltype(auto) getDynamicSizes() const
    {
      return subExpr.getDynamicSizes();
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    bool canAssign()
    {
      return false;
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=false;
    
    /////////////////////////////////////////////////////////////////
    
    //// Returns a component reducer on a different expression
    template <typename T>
    INLINE_FUNCTION
    decltype(auto) recreateFromExprs(T&& t) const
    {
      return compReduce<Combiner>(std::forward<T>(t));
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return compReduce<Rc,Combiner>(subExpr.getRef());	\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      CompReducer<Rc,
		  SameRefAs<_E,typename std::decay_t<_E>::template ReinterpretFund<NFund>>,
		  CompsList<C...>,
		  _Fund,
		  Combiner>;
    
    /////////////////////////////////////////////////////////////////
    
    /// Evaluate
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    Fund eval(const DerivedFromComp auto&...nTCs) const
    {
      /// Result
      Fund res;
      Combiner::setToInitialValue(res);
      
      compsLoop<CompsList<ReduceComp>>([this,&res,&nTCs...](const ReduceComp& rc) INLINE_ATTRIBUTE
      {
	/// First argument
	Combiner::reduce(res,this->subExpr(nTCs...,rc));
      },getDynamicSizes());
      
      return
	res;
    }
    
    /// Construct
    template <typename T>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    CompReducer(T&& arg)
      requires(std::is_same_v<std::decay_t<T>,std::decay_t<_E>>)
      : subExpr(std::forward<T>(arg))
    {
    }
    
    /// Copy constructor
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    CompReducer(const CompReducer& oth) :
      subExpr(oth.subExpr)
    {
    }
  };
  
  /// Reduces an expression over a component
  template <DerivedFromComp Rc,
	    typename Combiner,
	    DerivedFromNode _E>
  HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
  decltype(auto) compReduce(_E&& e)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    using Comps=
      TupleFilterAllTypes<typename E::Comps,CompsList<Rc>>;
    
    using Fund=
      typename E::Fund;
    
    return
      CompReducer<Rc,decltype(e),Comps,Fund,Combiner>(std::forward<_E>(e));
  }
  
#define PROVIDE_COMP_REDUCE_FUNCTOR(NAME,OP,INIT)		\
  namespace impl						\
  {								\
    /* Implements the reduction of a comp over operator NAME */	\
    struct _CompReduce## NAME ## Functor			\
    {								\
      template <typename Fund>					\
      INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB		\
      static void setToInitialValue(Fund& f)			\
      {								\
	f=INIT;							\
      }								\
      								\
      template <typename Fund>					\
      INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB		\
      static void reduce(Fund& out,				\
			 const Fund& in)			\
      {								\
	out OP##=in;						\
      }								\
    };								\
  }								\
								\
  /* Reduces the component C of expression E over operator NAME */\
  template <DerivedFromComp C,					\
	    DerivedFromNode E>					\
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
  decltype(auto) comp ## NAME (E&& e)				\
  {								\
    return compReduce<C,impl::_CompReduce ## NAME ##		\
		      Functor>(std::forward<E>(e));		\
  }
  
  PROVIDE_COMP_REDUCE_FUNCTOR(Sum,+,0);
  PROVIDE_COMP_REDUCE_FUNCTOR(Prod,*,1);
  
#undef PROVIDE_COMP_REDUCE_FUNCTOR
}

#endif
