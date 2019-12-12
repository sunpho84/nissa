#ifndef _STORAGE_HPP
#define _STORAGE_HPP

#include <base/memory_manager.hpp>
#include <base/metaprogramming.hpp>

namespace nissa
{
  /// Class to store the data
  template <typename Fund,   // Fundamental type
	    Size StaticSize, // Non-dynamic size
	    bool AllStatic>  // Store whether all components have static size
  struct TensStorage
  {
    /// Structure to hold dynamically allocated data
    struct DynamicStorage
    {
      /// Storage
      Fund* data;
      
      /// Construct allocating data
      DynamicStorage(const Size& dynSize)
      {
	data=memory_manager->provideAligned<Fund>(dynSize);
      }
      
      /// Destructor deallocating the memory
      ~DynamicStorage()
      {
	memory_manager->release(data);
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
    static constexpr bool stackAllocated=AllStatic and StaticSize*sizeof(Fund)<=MAX_STACK_SIZE;
    
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
