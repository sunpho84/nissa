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
  
  /// Tensor
  template <typename T,
	    typename...TC,
	    typename F>
  struct BaseTensor<T,TensorComps<TC...>,F> : Crtp<T>
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
    
    /// Full list of indices passed, not necessarily in the same order
    template <typename...TD>
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION
    const Fund& operator()(const TD&...td) const
    {
      return this->crtp().eval(td...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(operator(),CUDA_HOST_DEVICE);
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
