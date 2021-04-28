#ifndef _LOOKUPTABLE_HPP
#define _LOOKUPTABLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file lookupTable.hpp
///
/// \brief Lookup table that can be accessed both from host and device

/// Temporarily we nasty access and nasty allocate, then we well
/// register and add a method to initialize

#include <tensor/baseTensor.hpp>

namespace nissa
{
    template <typename TC,
	    typename Fund>
  struct LookupTable;

#define LOOKUP_TABLE LookupTable<TensorComps<TC...>,F>
#define BASE_TENSOR BaseTensor<LOOKUP_TABLE,TensorComps<TC...>,F>
  
  /// Tensor
  template <typename...TC,
	    typename F>
  struct LOOKUP_TABLE :
    BASE_TENSOR
  {
    using Base=BASE_TENSOR;
    
#undef BASE_TENSOR
#undef LOOKUP_TABLE
    
    /// Stored size
    Size storageSize;
    
    /// Storage allocated on host
    F* _host_storage;
    
#ifdef USE_CUDA
    /// Storage allocated on device
    F* _device_storage;
#endif
    
    /// Evaluate
    CUDA_HOST_DEVICE INLINE_FUNCTION const F& eval(const TC&...tc) const
    {
      return
#ifdef COMPILING_FOR_DEVICE
      _device_storage
#else
      _host_storage
#endif
	[this->indexComputer(tc...)];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(eval,CUDA_HOST_DEVICE);
    
    /// Mark if allocated
    bool allocated;
    
    /// Allocate the storage
    template <typename...TD>
    constexpr
    void allocate(const TensorComps<TD...>& td)
    {
      if(allocated)
	crash("Already allocated");
      
      this->indexComputer.setDynamicSizes(td);
      
      storageSize=this->indexComputer.maxVal();
      
      _host_storage=memoryManager<StorLoc::ON_CPU>()->template provide<F>(storageSize);
#ifdef USE_CUDA
      _device_storage=memoryManager<StorLoc::ON_GPU>()->template provide<F>(storageSize);
#endif
    }
    
    /// Allocate the storage when sizes are passed as a list of TensorComp
    template <typename...TDfeat>
    constexpr
    void allocate(const TensorCompFeat<TDfeat>&...tdFeat)
    {
      allocate(std::make_tuple(tdFeat.deFeat()...));
    }
    
    /// Destructor
    void dealloc()
    {
      if(allocated)
	{
	  memoryManager<StorLoc::ON_CPU>()->release(_host_storage);
#ifdef USE_CUDA
	  memoryManager<StorLoc::ON_GPU>()->release(_device_storage);
#endif
	}
      else
	crash("Not allocated");
      
      allocated=false;
      storageSize=0;
    }
    
    /// Copy on the device, after having been initilized
    void consolidate()
    {
#ifdef USE_CUDA
    cudaMemcpy(_device_storage,_host_storage,sizeof(F)*storageSize,cudaMemcpyHostToDevice);
#endif
    }
  };
}

#endif
