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
    explicit Tensor()
    {
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
#ifdef COMPILING_FOR_DEVICE
      if constexpr(SL==StorLoc::ON_CPU)
	__trap();
#else
      if constexpr(SL==StorLoc::ON_GPU)
	crash("Cannot access device memory from host");
#endif
      // /// Check that we are not accessing device memory from the host
      // constexpr bool accessDeviceMemoryOnHost=(SL==StorLoc::ON_GPU) and not CompilingForDevice;
      
      // /// Check that we are not accessing host memory on device
      // constexpr bool accessHostMemoryOnDevice=(SL==StorLoc::ON_CPU) and CompilingForDevice;
      
      // static_assert(not accessDeviceMemoryOnHost,"Cannot access device memory from host");
      // static_assert(not accessHostMemoryOnDevice,"Cannot access host memory from device");
      
      return storage[this->indexComputer(tc...)];
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
}

#endif
