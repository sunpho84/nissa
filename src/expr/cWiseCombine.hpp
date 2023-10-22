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
#include <expr/subNodes.hpp>
#include <metaprogramming/arithmeticTraits.hpp>
#include <tuples/tupleCat.hpp>
#include <tuples/uniqueTupleFromTuple.hpp>

namespace nissa
{
  PROVIDE_DETECTABLE_AS(CWiseCombiner);
  
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
    DynamicCompsProvider<CompsList<C...>>,
    DetectableAsCWiseCombiner,
    SubNodes<_E...>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    //static_assert(sizeof...(_E)==2,"Expecting 2 addends");
    
    IMPORT_SUBNODE_TYPES;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    // /// Execution space
    // static constexpr ExecSpace execSpace=
    //   commonExecSpace<std::remove_reference_t<_E>::execSpace...>();
    
    // static_assert(execSpace!=ExecSpace::HOST_DEVICE,"Cannot define coefficient wise combination in undefined exec space");
    
    /// List of dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<Comps>::DynamicComps;
    
    /// Sizes of the dynamic components
    const DynamicComps dynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
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
    
    /// Determine if coefficient-wise combination can be simdified - to be extended
    static constexpr bool simdifyCase()
    {
      return
	(SubNode<Is>::canSimdify and...) and
	std::is_same_v<typename SubNode<Is>::SimdifyingComp...>;
    }
    
    /// States whether the tensor can be simdified
    static constexpr bool canSimdify=
      simdifyCase();
    
    /// Type of the combining function
    using Comb=_Comb;
    
    // /// \todo improve
    
//     /// Components on which simdifying
//     using SimdifyingComp=typename SubNode<0>::SimdifyingComp;
    
// #define PROVIDE_SIMDIFY(ATTRIB)					\
//     /*! Returns a ATTRIB simdified view */			\
//     INLINE_FUNCTION						\
//     auto simdify() ATTRIB					\
//     {								\
//       return							\
// 	Comb::compute(SUBNODE(Is).simdify()...);		\
//     }
    
//     PROVIDE_SIMDIFY(const);
    
//     PROVIDE_SIMDIFY(/* non const */);
    
// #undef PROVIDE_SIMDIFY
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_RECREATE_FROM_EXPR(ATTRIB)			\
    /*! Returns a ATTRIB similar version */			\
    template <typename...T>					\
    INLINE_FUNCTION						\
    decltype(auto) recreateFromExprs(T&&...t) ATTRIB		\
    {								\
      return Comb::compute(std::forward<T>(t)...);		\
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
      return							\
	Comb::compute(SUBNODE(Is).getRef()...);			\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /// Gets the components for the I-th addend
    template <int I,
	      typename...Cs>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    static auto getCompsForAddend(const Cs&...cs)
    {
      return tupleGetSubset<typename SubNode<I>::Comps>(std::make_tuple(cs...));
    }
    
    /// Evaluate
    template <typename...Cs>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const Cs&...cs) const
    {
      return
	Comb::compute(std::apply(SUBNODE(Is),getCompsForAddend<Is>(cs...))...);
    }
    
    /// Construct
    template <typename...T>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    CWiseCombiner(const DynamicComps& dynamicSizes,
		   UNIVERSAL_CONSTRUCTOR_IDENTIFIER,
		   T&&...addends) :
      SubNodes<_E...>(addends...),
      dynamicSizes(dynamicSizes)
    {
    }
  };
  
  template <typename Comb,
	    typename..._E,
	    ENABLE_THIS_TEMPLATE_IF(isNode<_E> and...)>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  auto cWiseCombine(_E&&...e)
  {
    /// Computes the result components
    using Comps=
      UniqueTupleFromTuple<TupleCat<typename std::decay_t<_E>::Comps...>>;
    
    /// Determine the fundamental type of the result
    using Fund=
      decltype(Comb::compute(typename std::decay_t<_E>::Fund{}...));
    
    /// Resulting type
    using Res=
      CWiseCombiner<std::tuple<decltype(e)...>,
		    Comps,
		    Fund,
		    Comb>;
    
    /// Resulting dynamic components
    const auto dc=
      dynamicCompsCombiner<typename Res::DynamicComps>(e.getDynamicSizes()...);
    
    return
      Res(dc,UNIVERSAL_CONSTRUCTOR_CALL,std::forward<_E>(e)...);
  }
  
#define CATCH_BINARY_OPERATOR(OP,NAMED_OP)					\
  struct NAMED_OP							\
  {									\
    template <typename...Args>						\
    constexpr INLINE_FUNCTION						\
    static auto CUDA_HOST_AND_DEVICE compute(Args&&...s)		\
    {	/* forwarding is needed to preserve value category */		\
      return (std::forward<Args>(s) OP ...);				\
    }									\
  };  									\
  									\
  /*! Catch the OP operator */						\
  template <typename E1,						\
	    typename E2,						\
	    ENABLE_THIS_TEMPLATE_IF(isNode<E1> and isNode<E2>)>		\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
  auto operator OP(E1&& e1,						\
		   E2&& e2)						\
  {									\
									\
    return								\
      cWiseCombine<NAMED_OP>(std::forward<E1>(e1),std::forward<E2>(e2)); \
  }
  
  CATCH_BINARY_OPERATOR(+,plus);
  
  CATCH_BINARY_OPERATOR(-,minus);
  
  CATCH_BINARY_OPERATOR(/,divide);
  
  CATCH_BINARY_OPERATOR(%,modulo);
  
  CATCH_BINARY_OPERATOR(==,compare);
  
  CATCH_BINARY_OPERATOR(!=,differ);
  
  CATCH_BINARY_OPERATOR(and,andOp);
  
  CATCH_BINARY_OPERATOR(or,orOp);
  
  CATCH_BINARY_OPERATOR(xor,xorOp);
  
#undef CATCH_BINARY_OPERATOR
  
  /////////////////////////////////////////////////////////////////
  
#define CATCH_UNARY_OPERATOR(OP,NAMED_OP)				\
  struct NAMED_OP							\
  {									\
    template <typename Arg>						\
    constexpr INLINE_FUNCTION						\
    static auto CUDA_HOST_AND_DEVICE compute(Arg&& s)			\
    {	/* forwarding is needed to preserve value category */		\
      return OP std::forward<Arg>(s);					\
    }									\
  };  									\
  									\
  /*! Catch the OP operator */						\
  template <typename E,							\
	    ENABLE_THIS_TEMPLATE_IF(isNode<E>)>				\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
  auto operator OP(E&& e)						\
  {									\
									\
    return								\
      cWiseCombine<NAMED_OP>(std::forward<E>(e)); \
  }
  
  CATCH_UNARY_OPERATOR(+,uPlus);
  
  CATCH_UNARY_OPERATOR(-,uMinus);
  
  #undef CATCH_UNARY_OPERATOR
}

#endif
