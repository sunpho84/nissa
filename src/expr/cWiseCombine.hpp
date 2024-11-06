#ifndef _CWISE_COMBINE_HPP
#define _CWISE_COMBINE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/cWiseCombine.hpp

/// Combines multiple nodes with a specified function

#include <expr/comps.hpp>
#include <expr/conj.hpp>
#include <expr/node.hpp>
#include <metaprogramming/arithmeticTraits.hpp>
#include <tuples/tuple.hpp>
#include <tuples/tupleApply.hpp>
#include <tuples/tupleCat.hpp>
#include <tuples/uniqueTupleFromTuple.hpp>

namespace nissa
{
  PROVIDE_FEATURE(CWiseCombiner);
  
  /// CWiseCombiner
  ///
  /// Forward declaration to capture the components
  template <typename _E,
	    typename _Comps,
	    typename _Fund,
	    typename _Comb,
	    typename _Is=std::make_integer_sequence<int,std::tuple_size_v<_E>>>
  struct CWiseCombiner;
  
#define THIS								\
  CWiseCombiner<std::tuple<_E...>,CompsList<C...>,_Fund,_Comb,std::integer_sequence<int,Is...>>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// CWiseCombiner
  ///
  template <typename..._E,
	    typename...C,
	    typename _Fund,
	    typename _Comb,
	    int...Is>
  struct THIS :
    CWiseCombinerFeat,
    ManySubExprs<THIS>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Subexpressions
    Tuple<NodeRefOrVal<_E>...> subExprs;
    
    /// Type of I-th element
    template <size_t I>
    using SubExpr=
      std::decay_t<decltype(subExprs.template get<I>())>;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    /// Execution space
    static constexpr ExecSpace execSpace=
      (std::remove_reference_t<_E>::execSpace*...);
    
    /// List of dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<Comps>::DynamicComps;
    
    /// Sizes of the dynamic components
    const DynamicComps dynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    decltype(auto) getDynamicSizes() const
    {
      return dynamicSizes;
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    constexpr bool canAssign()
    {
      return false;
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=false;
    
    // /// Determine if coefficient-wise combination can be simdified - to be extended
    // static constexpr bool simdifyCase()
    // {
    //   return
    // 	(SubNode<Is>::canSimdify and...) and
    // 	std::is_same_v<typename SubNode<Is>::SimdifyingComp...>;
    // }
    
    // /// States whether the tensor can be simdified
    // static constexpr bool canSimdify=
    //   simdifyCase();
    
    /// Type of the combining function
    using Comb=_Comb;
    
    /////////////////////////////////////////////////////////////////
    
    /// Recreate from expressions
    template <typename... T>
    INLINE_FUNCTION
    auto recreateFromExprs(T&&...t) const
    {
      return Comb::compute(std::forward<T>(t)...);
    };
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return							\
	Comb::compute(subExprs.template get<Is>().getRef()...);	\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      CWiseCombiner<std::tuple<SameRefAs<_E,typename std::decay_t<_E>::template ReinterpretFund<NFund>>...>,
		    CompsList<C...>,
		    NFund,
		    _Comb,
		    std::integer_sequence<int,Is...>>;
    
    /// Gets the components for the I-th addend
    template <int I,
	      typename...Cs>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    static auto getCompsForSubexpr(const Cs&...cs)
    {
      return tupleGetSubset<typename SubExpr<I>::Comps>(std::make_tuple(cs...));
    }
    
    /// Evaluate
    template <typename...Cs>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    Fund eval(const Cs&...cs) const
    {
      return
	Comb::compute(tupleApply(subExprs.template get<Is>(),getCompsForSubexpr<Is>(cs...))...);
    }
    
    /// Construct
    template <DerivedFromNode...T>
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    CWiseCombiner(const DynamicComps& dynamicSizes,
		  const T&...addends) :
      subExprs{{addends}...},
      dynamicSizes(dynamicSizes)
    {
    }
    
    /// Copy constructor
    INLINE_FUNCTION constexpr
    CWiseCombiner(const CWiseCombiner& oth)=default;
    
    /// Move constructor
    INLINE_FUNCTION constexpr
    CWiseCombiner(CWiseCombiner&& oth)=default;
  };
  
  template <typename Comb,
	    DerivedFromNode..._E>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  auto cWiseCombine(const _E&...e)
  {
    /// Computes the result components
    using Comps=
      UniqueTupleFromTuple<TupleCat<typename std::decay_t<_E>::Comps...>>;
    
    /// Determine the fundamental type of the result
    using Fund=
      decltype(Comb::compute(typename std::decay_t<_E>::Fund{}...));
    
    /// Resulting type
    using Res=
      CWiseCombiner<std::tuple<const _E&...>,
		    Comps,
		    Fund,
		    Comb>;
    
    /// Resulting dynamic components
    const auto dc=
      dynamicCompsCombiner<typename Res::DynamicComps>(e.getDynamicSizes()...);
    
    return
      Res(dc,e...);
  }
  
#define CATCH_BINARY_OPERATOR(OP,NAMED_OP)				\
  namespace impl							\
  {									\
    struct _ ## NAMED_OP ## Functor					\
    {									\
      template <typename...Args>					\
      constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
      static auto compute(const Args&...s)				\
      {									\
	return (s OP ...);						\
      }									\
    };									\
  }									\
  									\
  /*! Catch the OP operator */						\
  template <DerivedFromNode E1,						\
	    DerivedFromNode E2>						\
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB				\
  auto operator OP(const E1& e1,					\
		   const E2& e2)					\
  {									\
									\
    return								\
      cWiseCombine<impl::_ ## NAMED_OP ## Functor>(e1,e2);		\
  }
  
  CATCH_BINARY_OPERATOR(+,Plus);
  
  CATCH_BINARY_OPERATOR(-,Minus);
  
  CATCH_BINARY_OPERATOR(/,Divide);
  
  CATCH_BINARY_OPERATOR(%,Modulo);
  
  CATCH_BINARY_OPERATOR(==,Compare);
  
  CATCH_BINARY_OPERATOR(!=,differ);
  
  CATCH_BINARY_OPERATOR(and,AndOp);
  
  CATCH_BINARY_OPERATOR(or,orOp);
  
  CATCH_BINARY_OPERATOR(xor,XorOp);
  
#undef CATCH_BINARY_OPERATOR
  
  /////////////////////////////////////////////////////////////////
  
#define CATCH_UNARY_OPERATOR(OP,NAMED_OP)				\
  namespace impl							\
  {									\
    struct _ ## NAMED_OP ## Functor					\
    {									\
      template <typename Arg>						\
      constexpr INLINE_FUNCTION						\
      static auto HOST_DEVICE_ATTRIB compute(Arg&& s)			\
      {									\
	return OP std::forward<Arg>(s);					\
      }									\
    };									\
  }									\
  									\
  /*! Catch the OP operator */						\
  template <DerivedFromNode E>						\
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB				\
  auto operator OP(const E& e)						\
  {									\
    return								\
      cWiseCombine<impl::_ ## NAMED_OP ## Functor>(e);			\
  }
  
  CATCH_UNARY_OPERATOR(+,uPlus);
  
  CATCH_UNARY_OPERATOR(-,uMinus);
  
#undef CATCH_UNARY_OPERATOR
}

#endif
