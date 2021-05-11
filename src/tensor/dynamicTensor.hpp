#ifndef _DYNAMIC_TENSOR_HPP
#define _DYNAMIC_TENSOR_HPP

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
#define TENSOR Tensor<TensorComps<TC...>,F,SL,TensorDynamicity::DYNAMIC_TENSOR>
#define BASE_TENSOR BaseTensor<TENSOR,TensorComps<TC...>,F>
  
  /// Tensor
  template <typename...TC,
	    typename F,
	    StorLoc SL>
  struct TENSOR :
    BASE_TENSOR
  {
    using Base=BASE_TENSOR;
    
#undef BASE_TENSOR
#undef TENSOR
    
    /// Fundamental type
    using Fund=typename Base::Fund;
    
    /// Components
    using Comps=typename Base::Comps;
    
    /// Get the I-th component
    template <int I>
    using Comp=
      std::tuple_element_t<I,Comps>;
    
    /// Storage
    Fund* storage;
    
    /// Storage size
    Size storageSize;
    
    bool allocated{false};
    
    /// Allocate the storage
    template <typename...TD>
    constexpr
    void allocate(const TensorComps<TD...>& td)
    {
      if(allocated)
	crash("Already allocated");
      
      this->indexComputer.setDynamicSizes(td);
      storageSize=this->indexComputer.maxVal();
      storage=memoryManager<SL>()->template provide<Fund>(storageSize);
    }
    
    /// Allocate the storage when sizes are passed as a list of TensorComp
    template <typename...TDfeat>
    constexpr
    void allocate(const TensorCompFeat<TDfeat>&...tdFeat)
    {
      allocate(std::make_tuple(tdFeat.deFeat()...));
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes
    template <typename...TD>
    constexpr
    explicit Tensor(const TensorComps<TD...>& td)
    {
      allocate(td);
    }
    
    /// Initialize the tensor without allocating
    constexpr
    Tensor()
    {
      allocated=false;
    }
    
    /// Initialize the tensor when sizes are passed as a list of TensorComp
    template <typename...TDfeat>
    constexpr
    explicit Tensor(const TensorCompFeat<TDfeat>&...td) :
      Tensor(std::make_tuple(td.deFeat()...))
    {
    }
    
    /// Destructor
    ~Tensor()
    {
      if(allocated)
	memoryManager<SL>()->release(storage);
      allocated=false;
    }
    
    // /// Move constructor
    // CUDA_HOST_DEVICE constexpr
    // Tensor(Tensor<TensorComps<TC...>,Fund,SL>&& oth) :
    //   dynamicSizes(oth.dynamicSizes),
    //   storage(std::move(oth.storage))
    // {
    // }
    
    // /// Move assignment
    // Tensor& operator=(Tensor&& oth)
    // {
    //   std::swap(dynamicSizes,oth.dynamicSizes);
    //   std::swap(storage,oth.data);
      
    //   return *this;
    // }
    
    // /// Copy constructor
    // Tensor(const Tensor& oth) :
    //   Tensor(oth.dynamicSizes)
    // {
    //   static_cast<Expr<Tensor,Comps>>(*this)=
    //   	static_cast<const Expr<Tensor,Comps>&>(oth);
    // }
    
    /// Check that we are accessing device vector only on device code
    void assertCorrectEvaluationStorage() const
    {
#ifdef COMPILING_FOR_DEVICE
      if constexpr(SL==StorLoc::ON_CPU)
	__trap();
#else
      if constexpr(SL==StorLoc::ON_GPU)
	crash("Cannot access device memory from host");
#endif
    }
    
#define PROVIDE_ORDERED_EVALUATOR(ATTRIB)					\
    /*! Evaluate, returning a reference to the fundamental type */	\
    template <typename...TD,						\
	      ENABLE_THIS_TEMPLATE_IF(std::is_same_v<TD,TC> && ...)>	\
    CUDA_HOST_DEVICE INLINE_FUNCTION					\
    ATTRIB Fund& orderedEval(const TD&...td) ATTRIB				\
    {									\
      assertCorrectEvaluationStorage();					\
									\
      return storage[this->indexComputer(td...)];			\
    }									\
    
    PROVIDE_ORDERED_EVALUATOR(/* not const */);
    
    PROVIDE_ORDERED_EVALUATOR(const);
    
#undef PROVIDE_ORDERED_EVALUATOR
  };
}

#endif
