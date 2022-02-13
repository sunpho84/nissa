#ifndef _SUM_HPP
#define _SUM_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file sum.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/nnaryExpr.hpp>

namespace nissa
{
#define THIS					\
  Summer<_Comps,_E1,_E2,_EvalTo>
  
#define NNEX					\
  NnaryExpr<THIS,				\
	    std::tuple<_E1,_E2>,		\
	    _Comps,				\
	    _EvalTo>
  
  DEFINE_FEATURE(Summer);
  
  /// Product of two expressions
  template <typename _Comps,
	    typename _E1,
	    typename _E2,
	    typename _EvalTo>
  struct Summer :
    SummerFeat<THIS>,
    NNEX
  {
    /// Import the nnary expression
    using NnEx=
      NNEX;
    
#undef UNEX
#undef THIS
    
    /// Components
    using Comps=
      _Comps;
    
    /// I-th Nested Expression type
    template <size_t I>
    using NestedExpr=
      std::decay_t<typename NnEx::template NestedExpr<I>>;
    
    /// Returned type when evaluating the expression
    using EvalTo=
      _EvalTo;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      EXPR_FLAG_MASKS::NONE;
    
    /// List of all dynamically allocated components
    using DynamicComps=
      GetDynamicCompsOfTensorComps<Comps>;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes;
    
    /// Dynamic sizes
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) getDynamicSizes() const
    {
      return
	dynamicSizes;
    }
    
    /// Evaluate one of the two addends, filtering the appropriate components
    ///
    /// Infer the components needed
    ///
    /// \todo since we do already a part of the eval at the expr stack
    /// level, we might call the orderedEval, but maybe we relook in
    /// the future
    template <int I,
	      typename TD,
	      typename...TC>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalAddend(const TD& td,
			      TensorComps<TC...>*) const
    {
      return
	this->template nestedExpr<I>()(std::get<TC>(td)...);
    }
    
    /// Evaluate one of the two addends, filtering the appropriate components
    ///
    /// Avoids relying on the current ability to disregard unneded
    /// components when evaluating an expression. If we keep the
    /// feature we might remove this call wrapper.
    ///
    /// \todo Should we move this functionality to the expr stack
    /// level?
    template <int I,
	      typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) evalAddend(const TD&...td) const
    {
      return
	_evalAddend<I>(std::make_tuple(td...),(typename NestedExpr<I>::Comps*)nullptr);
    }
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo eval(const TD&...td) const
    {
      return
	evalAddend<0>(td...)+
	evalAddend<1>(td...);
    }
    
    /// Construct
    template <typename E1,
	      typename E2,
	      typename...Td>
    Summer(E1&& e1,
	   E2&& e2,
	   const TensorComps<Td...>& dc) :
      NnEx(std::forward<E1>(e1),
	   std::forward<E2>(e2)),
      dynamicSizes(dc)
    {
    }
  };
  
  namespace internal
  {
    template <typename TC1,
	      typename TC2>
    struct _SumCompsComputer;
    
    template <typename...TC1,
	      typename...TC2>
    struct _SumCompsComputer<TensorComps<TC1...>,
			     TensorComps<TC2...>>
    {
      using type=
	UniqueTuple<TC1...,TC2...>;
    };
  }
  
  template <typename _E1,
	    typename _E2,
	    UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  auto sum(_E1&& e1,
	   _E2&& e2,
	   UNPRIORITIZE_DEFAULT_VERSION_ARGS)
  {
    UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK;
    
    /// Decayed type1
    using E1=
      std::decay_t<_E1>;
    
    /// Decayed type2
    using E2=
      std::decay_t<_E2>;
    
    /// First argument components
    using C1=
      typename E1::Comps;
    
    /// Second argument components
    using C2=
      typename E2::Comps;
    
    /// First expression fundamental type
    using EvalTo1=
      typename E1::EvalTo;
    
    /// Second expression fundamental type
    using EvalTo2=
      typename E2::EvalTo;
    
    /// Determine the fundamental type of the product
    using _EvalTo=
      decltype(EvalTo1()*EvalTo2());
    
    using Comps=
      typename internal::_SumCompsComputer<C1,C2>::type;
    
    /// Resulting type
    using Res=
      Summer<Comps,
	     decltype(e1),
	     decltype(e2),
	       _EvalTo>;
    
    /// Resulting dynamic components
    const auto& dc=
      dynamicCompsCombiner<typename Res::DynamicComps>(std::make_tuple(e1.getDynamicSizes(),e2.getDynamicSizes()));
    
    return
      Res(std::forward<_E1>(e1),
	  std::forward<_E2>(e2),
	  dc);
  }
  
  /// Catch the sum operator
  template <typename _E1,
	    typename _E2>
  auto operator+(const ExprFeat<_E1>& e1,
		 const ExprFeat<_E2>& e2)
  {
    return
      sum(e1.deFeat().crtp(),e2.deFeat().crtp());
  }
}

#endif
