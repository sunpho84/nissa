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
    
    /// Evaluate, returning a reference to the fundamental type
    ///
    /// Case in which the components are not yet correctly ordered
    template <typename...TD,
	      ENABLE_THIS_TEMPLATE_IF(sizeof...(TD)==sizeof...(TC))>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TD&...td) const
    {
      return eval(std::get<TC>(std::make_tuple(td...))...);
    }
    
    /// Evaluate, returning a reference to the fundamental type
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TC&...tc) const
    {
// #ifdef COMPILING_FOR_DEVICE
//       if constexpr(SL==StorLoc::ON_CPU)
// 	__trap();
// #else
//       if constexpr(SL==StorLoc::ON_GPU)
// 	crash("Cannot access device memory from host");
// #endif
      // /// Check that we are not accessing device memory from the host
      // constexpr bool accessDeviceMemoryOnHost=(SL==StorLoc::ON_GPU) and not CompilingForDevice;
      
      // /// Check that we are not accessing host memory on device
      // constexpr bool accessHostMemoryOnDevice=(SL==StorLoc::ON_CPU) and CompilingForDevice;
      
      // static_assert(not accessDeviceMemoryOnHost,"Cannot access device memory from host");
      // static_assert(not accessHostMemoryOnDevice,"Cannot access host memory from device");
      
      return this->crtp().storage[indexComputer(tc...)];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(eval,CUDA_HOST_DEVICE);
    
    /// Full list of indices passed, not necessarily in the same order
    template <typename...TD>
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION
    const Fund& operator()(const TD&...td) const
    {
      return eval(td...);
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
