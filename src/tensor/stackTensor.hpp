#ifndef _STACK_TENSOR_HPP
#define _STACK_TENSOR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tensor.hpp
///
/// \brief Implements all functionalities of tensors

#include <memory/memoryManager.hpp>
#include <memory/storLoc.hpp>
#include <tensor/baseTensor.hpp>

namespace nissa
{
#define TENSOR Tensor<TensorComps<TC...>,F,SL,TensorDynamicity::STACKED_TENSOR>
#define BASE_TENSOR BaseTensor<TENSOR,TensorComps<TC...>,F>
  
  /// Tensor
  template <typename...TC,
	    typename F,
	    StorLoc SL>
  struct TENSOR :
      BASE_TENSOR
  {
    // using Base=BASE_TENSOR;
    
#undef BASE_TENSOR
#undef TENSOR
    
    /// Fundamental type
    using Fund=F;
    
    /// Components
    using Comps=
      TensorComps<TC...>;
    
    /// Get the I-th component
    template <int I>
    using Comp=
      std::tuple_element_t<I,Comps>;
    
    /// Storage size
    static constexpr Size storageSize=
      IndexComputer<Comps>::staticPartMaxValue;
    
    /// Copy constructor
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    void nastyCopy(const Tensor& oth)
    {
      memcpy(this->storage,oth.storage,storageSize);
    }
    
    /// Actual storage
    Fund storage[storageSize];
    
    /// Evaluate, returning a reference to the fundamental type
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TC&...tc) const
    {
      return storage[this->indexComputer(tc...)];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(eval,CUDA_HOST_DEVICE);
  };
}

#endif
