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
    using Base=
      BASE_TENSOR;
    
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
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    void nastyCopy(const Tensor& oth)
    {
      memcpy(this->storage,oth.storage,storageSize);
    }
    
    /// Actual storage
    Fund storage[storageSize];
    
    /// THIS MIGHT BE MOVED TO A EVALUATOR FEATURE

#define PROVIDE_ORDERED_EVALUATOR(ATTRIB)				\
    /*! Evaluate, returning a reference to the fundamental type */	\
    template <typename...TD,						\
	      ENABLE_THIS_TEMPLATE_IF(std::is_same_v<TD,TC> && ...)>	\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION					\
    ATTRIB Fund& orderedEval(const TD&...td) ATTRIB			\
    {									\
      return storage[this->indexComputer(td...)];			\
    }									\
    
    PROVIDE_ORDERED_EVALUATOR(/* not const */);
    
    PROVIDE_ORDERED_EVALUATOR(const);
    
#undef PROVIDE_ORDERED_EVALUATOR
    
    using Base::operator=;
  };
}

#endif
