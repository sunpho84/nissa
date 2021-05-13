#ifndef _BASE_TENSOR_HPP
#define _BASE_TENSOR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file baseTensor.hpp
///
/// \brief Implements basic functionalities of tensors

#include <memory/memoryManager.hpp>
#include <memory/storLoc.hpp>
#include <metaProgramming/crtp.hpp>
#include <metaProgramming/refOrVal.hpp>
#include <tensor/expr.hpp>
#include <tensor/indexComputer.hpp>

namespace nissa
{
  /// Tensor with Comps components, of Fund fundamental type
  ///
  /// Forward definition to capture actual components
  template <typename T,
	    typename Comps,
	    typename Fund=double>
  struct BaseTensor;
  
  namespace details
  {
    constexpr ExprFlags _BaseTensorFlags=
      EXPR_FLAG_MASKS::EVAL_TO_REF|
      EXPR_FLAG_MASKS::STORE_BY_REF|
      EXPR_FLAG_MASKS::NEEDS_ORDERED_COMPS_TO_EVAL;
  }
  
  /// Tensor
  template <typename T,
	    typename...TC,
	    typename F>
  struct BaseTensor<T,TensorComps<TC...>,F>
    : Expr<T,TensorComps<TC...>,F,
	   details::_BaseTensorFlags> // Here we pass the external expression on purpose
  {
    /// Fundamental type
    using Fund=F;
    
    /// Components
    using Comps=
      TensorComps<TC...>;
    
    /// Get the I-th component
    template <int I>
    using Comp=
      std::tuple_element_t<I,Comps>;
    
    /// Index computer type
    using IC=IndexComputer<TensorComps<TC...>>;
    
    /// Index computer
    IC indexComputer;
    
    /// Dynamic sizes
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) getDynamicSizes() const
    {
      return
	indexComputer.dynamicSizes;
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  enum class TensorDynamicity{STACKED_TENSOR,DYNAMIC_TENSOR};
  
  namespace details
  {
    template <typename TC,
	      typename Fund,
	      StorLoc SL>
    struct _TensorStackedDecider
    {
      using IC=IndexComputer<TC>;
      
      static constexpr auto staticSize=
	IC::staticPartMaxValue*sizeof(Fund);
      
      /// Threshold beyond which allocate dynamically in any case
      static constexpr
      Size MAX_STACK_SIZE=2304;
      
      static constexpr bool isSmall=
	staticSize<MAX_STACK_SIZE;
      
      static constexpr bool hasNoDynamicComps=
	 IC::allCompsAreStatic;
      
      static constexpr bool compilingForSameLoc=
	     ((CompilingForDevice==true  and SL==StorLoc::ON_GPU) or
	      (CompilingForDevice==false and SL==StorLoc::ON_CPU));
      
      static constexpr
      TensorDynamicity tensorDynamicity=
	(isSmall and hasNoDynamicComps and compilingForSameLoc)?
	TensorDynamicity::STACKED_TENSOR:
	TensorDynamicity::DYNAMIC_TENSOR;
    };
  }
  
  template <typename TC,
	    typename Fund=double,
	    StorLoc SL=DefaultStorage,
	    TensorDynamicity Dynamicity=details::_TensorStackedDecider<TC,Fund,SL>::tensorDynamicity>
  struct Tensor;
}

#endif
