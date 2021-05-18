#ifndef _PROD_HPP
#define _PROD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file prod.hpp

#include <metaProgramming/universalReference.hpp>
#include <tensor/conj.hpp>
#include <tensor/loopOnAllComponents.hpp>
#include <tensor/nnaryExpr.hpp>

namespace nissa
{
  namespace internal
  {
    /// Classifies the components, determining which one are visible or contracted
    ///
    /// Internal implementation, forward declaration
    template <typename TT,
	      typename TP1,
	      typename TP2>
    struct _ProdCompsClassifierImpl;
    
    /// Classifies the components, determining which one are visible or contracted
    template <typename...S,
	      RwCl...RCA,
	      int...Which,
	      typename TP1,
	      typename TP2>
    struct _ProdCompsClassifierImpl<TensorComps<TensorComp<S,RCA,Which>...>,
				    TP1,
				    TP2>
    {
      /// Detect if component with signature _S and index _Which is contracted
      template <typename  _S,
		int _Which>
      static constexpr bool contract=
	tupleHasType<TP1,TensorComp<_S,CLN,_Which>,1> and
	tupleHasType<TP2,TensorComp<_S,ROW,_Which>,1>;
      
      /// Classifies the component TC
      ///
      /// Forward declaration
      template <typename TC>
      struct _Classify;
      
      /// Classifies a row component
      template <typename _S,
		int _Which>
      struct _Classify<TensorComp<_S,ROW,_Which>>
      {
	/// This component, in row format
	using ThisRow=
	  TensorComp<_S,ROW,_Which>;
	
	/// This component, in cln format
	using ThisCln=
	  TensorComp<_S,CLN,_Which>;
	
	/// Determine if the result has the row component
	static constexpr bool resHasRow=
	  tupleHasType<TP1,ThisRow,1> or (tupleHasType<TP2,ThisRow,1> and not contract<_S,_Which>);
	
	/// Determine if the result has the cln component
	static constexpr bool resHasCln=
	  tupleHasType<TP2,ThisCln,1> or (tupleHasType<TP1,ThisCln,1> and not contract<_S,_Which>);
	
	/// Put together the two components
	using VisibleComps=
	  TupleCat<std::conditional_t<resHasRow,TensorComps<ThisRow>,TensorComps<>>,
		   std::conditional_t<resHasCln,TensorComps<ThisCln>,TensorComps<>>>;
      };
      
      /// Classifies a non row-cln component
      template <typename _S,
		int _Which>
      struct _Classify<TensorComp<_S,ANY,_Which>>
      {
	
	using VisibleComps=
	  TensorComps<TensorComp<_S,ANY,_Which>>;
      };
      
      /// Cat all visible components
      using VisibleComps=
	TupleCat<typename _Classify<TensorComp<S,RCA,Which>>::VisibleComps...>;
      
      /// Cat all contracted components
      using ContractedComps=
	TupleCat<std::conditional_t<contract<S,Which>,TensorComps<TensorComp<S,ROW,Which>>,TensorComps<>>...>;
    };
    
    template <typename TC1,
	      typename TC2>
    using _ProdCompsComputer=
      _ProdCompsClassifierImpl<IndependentComponents<TC1,TC2>,
			   TC1,TC2>;
  }
  
  /////////////////////////////////////////////////////////////////
  
  DEFINE_FEATURE(Producer);
  
#define THIS					\
  Producer<_ContractedComps,_E1,_E2,_Comps,_EvalTo>
  
#define NNEX					\
  NnaryExpr<THIS,				\
	    std::tuple<_E1,_E2>,		\
	    _Comps,				\
	    _EvalTo>
  
  /// Product of two expressions
  template <typename _ContractedComps,
	    typename _E1,
	    typename _E2,
	    typename _Comps,
	    typename _EvalTo>
  struct Producer :
    ProducerFeat<THIS>,
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
    
    /// Detect complex component
    static constexpr bool isComplProd=
      tupleHasType<Comps,ComplId>;
    
    /// Evaluate the I-th argument
    template <size_t I,
	      typename...Args>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalIthArgFilteredList(const Args&...args) const
    {
      return
	this->template nestedExpr<I>()(args...);
    }
    
    /// Evaluate the I-th argument, expanding the tuple containing the arguments into a list
    template <size_t I,
	      typename...Args>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalIthArgFiltered(const TensorComps<Args...> arg) const
    {
      return
	_evalIthArgFilteredList<I>(std::get<Args>(arg)...);
    }
    
    /// Evaluate the I-th argument, gathering only the components that need actually be passed to the subexpr
    template <size_t I,
	      typename...Args>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalIthArg(const TensorComps<Args...> args) const
    {
      /// I-th Nested Expression type
      using NestedExprI=
	std::decay_t<typename NnEx::template NestedExpr<I>>;
      
      /// Components of the subexpression
      using SubComps=
	typename NestedExprI::Comps;
      
      return
	_evalIthArgFiltered<I>(tupleGetSubset<SubComps>(args));
    }
    
    /// Decide the strategy to evaluate
    struct _EvalStrategy
    {
      /// Possible strategies
      enum{NON_COMPL_PROD,COMPL_PROD};
      
      /// Decides the strategy for a given component
      using Get=
	std::integral_constant<int,
			       ((not isComplProd)?
			       NON_COMPL_PROD:
			       COMPL_PROD)>*;
      
      DECLARE_STRATEGY(NonComplProd,NON_COMPL_PROD);
      
      DECLARE_STRATEGY(ComplProd,COMPL_PROD);
    };
    
    /// Evaluate for non complex expressions
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo _eval(typename _EvalStrategy::NonComplProd,
		 const TD&...td) const
    {
      /// Result
      EvalTo res(0);
      
      loopOnAllComponents<_ContractedComps>(TensorComps<>(),
					    [this,&res,&td...](const auto&...args) INLINE_ATTRIBUTE
      {
	/// Put together all arguments
	const auto fullArgs=
	  std::make_tuple(td...,args...,args.transp()...);
	
	res+=
	  this->_evalIthArg<0>(fullArgs)*
	  this->_evalIthArg<1>(fullArgs);
	});
      
      return
	res;
    }
    
    /// Evaluate for complex expressions
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo _eval(typename _EvalStrategy::ComplProd,
		 const TD&...td) const
    {
      /// Result
      EvalTo res(0);
      
      loopOnAllComponents<_ContractedComps>(getDynamicSizes(),
					    [this,&res,&td...](const auto&...args) INLINE_ATTRIBUTE
      {
	/// Gets the complex index
	const ComplId cId=
	  std::get<ComplId>(std::make_tuple(td...));
	
	/// Put together all arguments
	auto fullArgs=
	  std::make_tuple(td...,args...,args.transp()...);
	
	std::get<ComplId>(fullArgs)=
	  Re;
	
	/// Compute the real part of argument 0
	const auto re0=
	  this->_evalIthArg<0>(fullArgs);
	  
	/// Compute the real part of argument 1
	const auto re1=
	  this->_evalIthArg<1>(fullArgs);
	
	std::get<ComplId>(fullArgs)=
	  Im;
	
	/// Compute the imaginary part of argument 0
	const auto im0=
	  this->_evalIthArg<0>(fullArgs);
	  
	/// Compute the imaginary part of argument 1
	const auto im1=
	  this->_evalIthArg<1>(fullArgs);
	
	if(cId==Re)
	  {
	    res+=
	      re0*re1;
	    res-=
	      im0*im1;
	  }
	else
	  {
	    res+=
	      re0*im1;
	    res+=
	      im0*re1;
	  }
	});
      
      return
	res;
    }
    
    /// Evaluate dispatching the correct strategy
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo eval(const TD&...td) const
    {
      return
	_eval(typename _EvalStrategy::Get(),
	      td...);
    }
    
    /// Construct
    template <typename E1,
	      typename E2,
	      typename...Td>
    Producer(E1&& e1,
	     E2&& e2,
	     const TensorComps<Td...>& dc) :
      NnEx(std::forward<E1>(e1),
	   std::forward<E2>(e2)),
      dynamicSizes(dc)
    {
    }
  };
  
  template <typename _E1,
	    typename _E2,
	    UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  auto prod(_E1&& e1,
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
    
    /// Computes the product components
    using PCC=
      internal::_ProdCompsComputer<C1,C2>;
    
    /// Gets the visible comps
    using VisibleComps=
      typename PCC::VisibleComps;
    
    /// Gets the contracted comps
    using ContractedComps=
      typename PCC::ContractedComps;
    
    /// First expression fundamental type
    using EvalTo1=
      typename E1::EvalTo;
    
    /// Second expression fundamental type
    using EvalTo2=
      typename E2::EvalTo;
    
    /// Determine the fundamental type of the product
    using _EvalTo=
      decltype(EvalTo1()*EvalTo2());
    
    /// Resulting type
    using Res=
      Producer<ContractedComps,
	       decltype(e1),
	       decltype(e2),
	       VisibleComps,
	       _EvalTo>;
    
    /// Resulting dynamic components
    const auto& dc=
      dynamicCompsCombiner<typename Res::DynamicComps>(std::make_tuple(e1.getDynamicSizes(),e2.getDynamicSizes()));
    
    return
      Res(std::forward<_E1>(e1),
	  std::forward<_E2>(e2),
	  dc);
  }
  
  /// Catch the product operator
  template <typename _E1,
	    typename _E2>
  auto operator*(const ExprFeat<_E1>& e1,
		 const ExprFeat<_E2>& e2)
  {
    return
      prod(e1.deFeat().crtp(),e2.deFeat().crtp());
  }
}

#endif
