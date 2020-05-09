#ifndef _STORAGE_HPP
#define _STORAGE_HPP

#include <base/memory_manager.hpp>
#include <base/metaprogramming.hpp>
#include <field/component.hpp>

namespace nissa
{
  /// Position where to store the data: device or host
  enum class TensStorageLocation{ON_CPU
#ifdef USE_CUDA
				 ,ON_GPU
#endif
  };
  
  /// Wraps the memory manager
  template <TensStorageLocation>
  struct MemoryManageWrapper;
  
  /// Use memory manager
  template <>
  struct MemoryManageWrapper<TensStorageLocation::ON_CPU>
  {
    static auto& get()
    {
      return cpu_memory_manager;
    }
  };
  
#ifdef USE_CUDA
  /// Use memory manager
  template <>
  struct MemoryManageWrapper<TensStorageLocation::ON_GPU>
  {
    static auto& get()
    {
      return gpu_memory_manager;
    }
  };
#endif
  
  /// Class to store the data
  template <typename Fund,           // Fundamental type
	    Size StaticSize,         // Size konwn at compile time
	    TensStorageLocation SL>  // Location where to store data
  struct TensStorage
  {
    /// Structure to hold dynamically allocated data
    struct DynamicStorage
    {
      /// Gets the appropriate memory manager
      static auto memory_manager()
      {
	return MemoryManageWrapper<SL>::get();
      }
      
      /// Storage
      Fund* data;
      
      /// Construct allocating data
      DynamicStorage(const Size& dynSize)
      {
	memory_manager()->template provide<Fund>(dynSize);
      }
      
      /// Destructor deallocating the memory
      ~DynamicStorage()
      {
	memory_manager()->release(data);
      }
    };
    
    /// Structure to hold statically allocated data
    struct StackStorage
    {
      /// Storage
      Fund data[StaticSize];
      
      /// Constructor: since the data is statically allocated, we need to do nothing
      StackStorage(const Size&)
      {
      }
    };
    
    /// Threshold beyond which allocate dynamically in any case
    static constexpr Size MAX_STACK_SIZE=2304;
    
    /// Decide whether to allocate on the stack or dynamically
    static constexpr bool stackAllocated=StaticSize!=DYNAMIC and StaticSize*sizeof(Fund)<=MAX_STACK_SIZE;
    
    /// Actual storage class
    using ActualStorage=std::conditional_t<stackAllocated,StackStorage,DynamicStorage>;
    
    /// Storage of data
    ActualStorage data;
    
    /// Forbids copy
    TensStorage(const TensStorage&) =delete;
    
    /// Construct taking the size to allocate
    TensStorage(const Size& size) ///< Size to allocate
      : data(size)
    {
    }
    
    /// Single component access via subscribe operator
    template <typename T>                       // Subscribed component type
    const auto& operator[](const T& t) const  ///< Subscribed component
    {
      return data.data[t];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(operator[]);
  };
}

#endif
