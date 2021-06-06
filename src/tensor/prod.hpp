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
  
  namespace internal
  {
    /// Gets the mapping between product components and one of its factor
    ///
    /// Forward implementation
    template <typename FCs, // Factor components
	      typename ECs> // Excluded (contracted) components
    struct _CompsRemappingForProductFactor;
    
    /// Gets the mapping between product components and one of its factor
    ///
    /// Invariant are passed as structure template arguments, while
    /// the product arguments are scanned from the Getter substruct
    template <typename...FactorComps,
	      typename...ExcludedComps>
    struct _CompsRemappingForProductFactor<TensorComps<FactorComps...>,
					   TensorComps<ExcludedComps...>>
    {
      /// Tag to be used for "excluded" is the same of not found, that is, val beyond the end
      static constexpr size_t notFoundTag=
	sizeof...(FactorComps);
      
      /// Check if a given component is in the excluded list
      template <typename PComp>
      static constexpr bool isExcluded=
	(std::is_same_v<ExcludedComps,PComp>||...);
      
      /// Take the position of the given component in the factor components
      template <typename PComp>
      static constexpr size_t pos=
	firstOccurrenceOfTypeInList<PComp,FactorComps...>;
      
      /// Returns the position, or value past end if not found or excluded
      template <typename PComp>
      static constexpr size_t value=
	isExcluded<PComp>?
	notFoundTag:
	pos<PComp>;
      
      /// Produce the resulting list
      ///
      /// Forward declaration
      template <typename PComps>
      struct _Getter;
      
      /// Produce the resulting list
      template <typename...PComps>
      struct _Getter<TensorComps<PComps...>>
      {
	/// Resulting type
	using type=
	  std::index_sequence<value<PComps>...>;
      };
    };
    
    /// Gets the mapping between product components and its factors
    ///
    /// Internal implementation, forward declaration
    template <typename PCs,     // Product components
	      typename Fs,      // Factors
	      typename CCs,     // Excluded (contracted) components
	      bool isComplProd, // Takes note if this is complex product
	      typename Is=      // Integer sequence to scan both args
	      std::index_sequence<0,1>>
    struct _GetCompsMeldBarriersForProduct;
    
    /// Gets the mapping between product components and its factors
    ///
    /// Internal implementation
    /// The original implementation, with unrolled factors, can be
    /// found on testNewExpr.cpp at commit 481b4908
    template <typename PCs,     // Product components
	      typename Fs,      // Factors
	      typename CCs,     // Excluded (contracted) components
	      bool isComplProd, // Takes note if this is complex product
	      size_t...Is>      // Integer sequence to scan both args
    struct _GetCompsMeldBarriersForProduct<PCs,Fs,CCs,isComplProd,std::index_sequence<Is...>>
    {
      /// I-th Factor
      template <size_t I>
      using F=
	typename std::tuple_element_t<I,Fs>;
      
      /// I-th Factor Components
      template <size_t I>
      using C=
	typename F<I>::Comps;
      
      /// I-th Factor meld barriers
      template <size_t I>
      using M=
	typename F<I>::CompsMeldBarriers;
      
      /// Components to exclude
      ///
      /// N.B: Second arguments needs to exclude the contracted
      /// components, first factor needs to exclude transposed ones
      template <size_t I>
      using EC=
	std::conditional_t<I==0,
			   TransposeMatrixTensorComps<CCs>,
			   CCs>;
      
      /// I-th remapping
      template <size_t I>
      using Remap=
	typename _CompsRemappingForProductFactor<C<I>,EC<I>>::template _Getter<PCs>::type;
      
      /// Resulting barriers from remapping
      using RemapBarr=
	typename internal::_GetTensorCompsMeldBarriersFromCompsRemapping<PCs,
									 std::tuple<C<Is>...>,
									 std::tuple<M<Is>...>,
									 std::tuple<Remap<Is>...>>::type;
      
      /// Position of the complex component
      static constexpr size_t complPos=
	firstOccurrenceOfTypeInTuple<ComplId,PCs>;
      
      /// Insert a barrier to avoid melding complex component
      static constexpr auto _insertComplBarr(std::bool_constant<true>)
      {
	using ComplProdBarr=
	  std::index_sequence<complPos,complPos+1>;
	
	using type=
	  CompsMeldBarriersInsert<RemapBarr,ComplProdBarr>;
	
	return
	  type{};
      }
      
      /// Do not insert a barrier to avoid melding complex component
      static constexpr auto _insertComplBarr(std::bool_constant<false>)
      {
	using type=
	  RemapBarr;
	
	return type{};
      };
      
      /// Iclude a barrier to avoid melding complex component if this is a complex product
      using type=
	decltype(_insertComplBarr(std::bool_constant<isComplProd>()));
   };
  }
  
  /// Gets the mapping between product components and its factors
  template <typename PCs,     // Product components
	    typename F1,      // First factor
	    typename F2,      // Second factor
	    typename CCs,     // Excluded (contracted) components
	    bool isComplProd> // Takes note if this is complex product
  using GetCompsMeldBarriersForProduct=
    typename internal::_GetCompsMeldBarriersForProduct<PCs,std::tuple<F1,F2>,CCs,isComplProd>::type;
  
  /////////////////////////////////////////////////////////////////
  
  DEFINE_FEATURE(Producer);
  
#define THIS					\
  Producer<_ContractedComps,_E1,_E2,_Comps,_CompsMeldBarriers,_isComplProd,_EvalTo>
  
#define NNEX					\
  NnaryExpr<THIS,				\
	    std::tuple<_E1,_E2>,		\
	    _Comps,				\
	    _CompsMeldBarriers,		\
	    _EvalTo>
  
  /// Product of two expressions
  template <typename _ContractedComps,
	    typename _E1,
	    typename _E2,
	    typename _Comps,
	    typename _CompsMeldBarriers,
	    bool _isComplProd,
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
    
    /// I-th Nested Expression type
    template <size_t I>
    using NestedExpr=
      std::decay_t<typename NnEx::template NestedExpr<I>>;
    
    template <size_t I>
    using NonContractedSubComps=
      TupleFilterAllTypes<typename NestedExpr<I>::Comps,
			  std::conditional_t<I==0,TransposeMatrixTensorComps<_ContractedComps>,_ContractedComps>>;
    
    /// Returned type when evaluating the expression
    using EvalTo=
      _EvalTo;
    
    /// List of all contracted components
    using ContractedComps=
      _ContractedComps;
    
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
    
    /// Takes note wheter this is complex product
    static constexpr bool isComplProd=
      _isComplProd;
    
    /// Evaluate the I-th argument, expanding the tuple containing the arguments into a list
    template <typename...Comps,
	      typename Arg>
    static CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalArgFiltered(const TensorComps<Comps...> comps,
				       Arg&& arg)
    {
      return
	arg(std::get<Comps>(comps)...);
    }
    
    /// Evaluate the I-th argument, gathering only the components that need actually be passed to the subexpr
    template <typename NonContractedCompsOfArg,
	      typename NonContractedComps,
	      typename ContractedCompsToPassToArg,
	      typename Arg>
    static CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalArgAfterFiltering(const NonContractedComps& nonContractedComps,
			    const ContractedCompsToPassToArg contractedCompsToPassToArg,
			    Arg&& arg)
    {
      return
	_evalArgFiltered(std::tuple_cat(tupleGetSubset<NonContractedCompsOfArg>(nonContractedComps),contractedCompsToPassToArg),std::forward<Arg>(arg));
    }
    
    /// Returns the contracted components in the format needed to evaluate the argument
    ///
    /// If the argument is the first one, we need to transpose the
    /// components, if it is the second, we need not
    template <size_t I,
	      typename...ContractedComps>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    static auto getContractedCompsForArg(const ContractedComps...contractedComps)
    {
      if constexpr(I==0)
	return
	  std::make_tuple(contractedComps.transp()...);
      else
	return
	  std::make_tuple(contractedComps...);
    }
    
    /// Evaluate the required argument, getting the contracted components in the proprer shape
    template <size_t I,
	      typename NonContractedComps,
	      typename...ContractedComps>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) _evalArg(const NonContractedComps& nonContractedComps,
			    const ContractedComps...contractedComps) const
    {
      return
	_evalArgAfterFiltering<NonContractedSubComps<I>>(nonContractedComps,
							 getContractedCompsForArg<I>(contractedComps...),
							 this->template nestedExpr<I>());
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
      
      DECLARE_DISPATCH_STRATEGY(NonComplProd,NON_COMPL_PROD);
      
      DECLARE_DISPATCH_STRATEGY(ComplProd,COMPL_PROD);
    };
    
    DECLARE_COMPONENT(SubExpr,int,2);
    
    /// Evaluate for non complex expressions
    template <typename...NCCs>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo _eval(typename _EvalStrategy::NonComplProd,
		 const NCCs&...nCCs) const
    {
      /// Result
      EvalTo res(0);
      
      loopOnAllComponents<_ContractedComps>(TensorComps<>(),
					    [this,&res,nCCs=std::make_tuple(nCCs...)](const auto&...cCs) INLINE_ATTRIBUTE
      {
	/// First argument
	const auto f=
	  this->_evalArg<0>(nCCs,cCs...);
	
	/// Second argument
	const auto s=
	  this->_evalArg<1>(nCCs,cCs...);
	
	res+=
	  f*s;
	});
      
      return
	res;
    }
    
    /// Replace the complex index in the list with the passed one
    ///
    /// Internal implementation
    template <size_t...HeadPos,    // Positions before complex id
	      size_t ComplPos,     // Position of complex id
	      size_t...TailPos,    // Positions after complex id
	      typename NCCs>       // Non contracted comps
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    auto _getNonContractedCompsForComplProdPart(std::index_sequence<HeadPos...>,
						std::index_sequence<ComplPos>,
						std::index_sequence<TailPos...>,
						const ComplId& cId,
						const NCCs& nCCs) const
    {
      return
	std::make_tuple(std::get<HeadPos>(nCCs)...,cId,std::get<TailPos+ComplPos+1>(nCCs)...);
    }
    
    /// Replace the complex index in the list with the passed one
    template <typename NCCs>      // Non contracted comps
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    auto _getNonContractedCompsForComplProdPart(const ComplId& cId,
						const NCCs& nCCs) const
    {
      /// Get the position of complex index
      constexpr size_t complPos=
	firstOccurrenceOfTypeInTuple<ComplId,NCCs>;
	
      return
	_getNonContractedCompsForComplProdPart(std::make_index_sequence<complPos>{},
					       std::index_sequence<complPos>{},
					       std::make_index_sequence<std::tuple_size_v<NCCs>-1-complPos>{},
					       cId,
					       nCCs);
    }
    
    /// Evaluate for complex expressions, dispatching the components to real and imaginary part
    template <typename NCCs,     // Non contracted comps
	      typename...CCs>      // Contracted comps
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo _evalComplProd(const NCCs& nCCs,
			  const CCs&...cCs) const
    {
      /// Non contracted components needed to evaluate the real part
      const auto nCCsRe=
	_getNonContractedCompsForComplProdPart(Re,nCCs);
      
      /// Non contracted components needed to evaluate the imag part
      const auto nCCsIm=
	_getNonContractedCompsForComplProdPart(Im,nCCs);
      
      /// Compute the real part of argument 0
      decltype(auto) re0=
	_evalArg<0>(nCCsRe,cCs...);
      
      /// Compute the real part of argument 1
      decltype(auto) re1=
	_evalArg<1>(nCCsRe,cCs...);
      
      /// Compute the imag part of argument 0
      decltype(auto) im0=
	_evalArg<0>(nCCsIm,cCs...);
      
      /// Compute the imag part of argument 1
      decltype(auto) im1=
	_evalArg<1>(nCCsIm,cCs...);
      
      /// Get the position of complex index
      const ComplId& cId=
	std::get<ComplId>(nCCs);
      
      if(cId==Re)
	return
	  re0*re1
	  -
	  im0*im1;
      else
	return
	  re0*im1
	  +
	  im0*re1;
    }
    
    /// Evaluate for complex expressions
    template <typename...NCCs>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo _eval(typename _EvalStrategy::ComplProd,
		 const NCCs&...nCCs) const
    {
      /// Result
      EvalTo res(0);
      
      loopOnAllComponents<_ContractedComps>(getDynamicSizes(),
					    [this,&res,&nCCs...](const auto&...cCs) INLINE_ATTRIBUTE
      {
	res+=
	  this->_evalComplProd(std::make_tuple(nCCs...),
			       cCs...);
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
    
    /// Detect complex product
    static constexpr bool isComplProd=
      tupleHasType<C1,ComplId> and
      tupleHasType<C2,ComplId>;
    
    /// Fusable comps
    using CompsMeldBarriers=
      GetCompsMeldBarriersForProduct<VisibleComps,E1,E2,ContractedComps,isComplProd>;
    
    /// Resulting type
    using Res=
      Producer<ContractedComps,
	       decltype(e1),
	       decltype(e2),
	       VisibleComps,
	       CompsMeldBarriers,
	       isComplProd,
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
