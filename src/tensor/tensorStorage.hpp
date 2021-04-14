#ifndef _TENSOR_STORAGE_HPP
#define _TENSOR_STORAGE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <memory/memoryManager.hpp>
#include <memory/storLoc.hpp>
#include <tensor/component.hpp>

/// \file tensorStorage.hpp

namespace nissa
{
  /// Stackability
  enum class Stackable{CANNOT_GO_ON_STACK,MIGHT_GO_ON_STACK};
  
  /// Basic storage, to use to detect storage
  template <typename T>
  struct BaseTensorStorage
  {
    PROVIDE_DEFEAT_METHOD(T);
  };
  
  /// Class to store the data
  template <typename Fund,            // Fundamental type
	    Size StaticSize,          // Size konwn at compile time
	    StorLoc SL,               // Location where to store data
	    Stackable IsStackable=    // Select if can go or not on stack
	    Stackable::MIGHT_GO_ON_STACK>
  struct TensorStorage
  {
    /// Structure to hold dynamically allocated data
    struct DynamicStorage
    {
      /// Hold info if it is a reference
      bool isRef;
      
      /// Storage
      Fund* innerStorage;
      
      /// Allocated size
      Size dynSize;
      
      /// Returns the size
      constexpr Size getSize()
	const
      {
	return
	  dynSize;
      }
      
      /// Returns the pointer to data
      CUDA_HOST_DEVICE
      Fund getDataPtr() const
      {
	return innerStorage;
      }
      
      PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(getDataPtr,CUDA_HOST_DEVICE);
      
      /// Allocate the given size
      void allocate(const Size& extSize)
      {
	if(not isRef and dynSize!=0)
	  crash("Already allocated, with memory size: %ld",dynSize);
	
	isRef=true;
	dynSize=extSize;
	innerStorage=memoryManager<SL>()->template provide<Fund>(dynSize);
      }
      
      /// Default constructor
      DynamicStorage() :
	isRef(true),
	innerStorage(nullptr),
	dynSize(0)
      {
      }
      
      /// Construct allocating data
      DynamicStorage(const Size& dynSize) :
	DynamicStorage()
      {
	allocate(dynSize);
      }
      
      /// "Copy" constructor, actually taking a reference
      CUDA_HOST_DEVICE
      DynamicStorage(const DynamicStorage& oth) :
	isRef(true),
	dynSize(oth.dynSize),
	innerStorage(oth.innerStorage)
      {
      }
      
      /// Create a reference starting from a pointer
      CUDA_HOST_DEVICE
      DynamicStorage(Fund* oth,
		     const Size& dynSize) :
	isRef(true),
	innerStorage(oth),
	dynSize(dynSize)
      {
      }
      
      /// Move constructor
      CUDA_HOST_DEVICE
      DynamicStorage(DynamicStorage&& oth) :
	isRef(oth.isRef),
	innerStorage(oth.innerStorage),
	dynSize(oth.dynSize)
      {
	oth.isRef=true;
	oth.innerStorage=nullptr;
	oth.dynSize=0;
      }
      
      /// Move assignment
      DynamicStorage& operator=(DynamicStorage&& oth)
      {
	std::swap(isRef,oth.isRef);
	std::swap(innerStorage,oth.innerStorage);
	std::swap(dynSize,oth.dynSize);
	
	return *this;
      }
      
#ifndef COMPILING_FOR_DEVICE
      /// Destructor deallocating the memory
      ~DynamicStorage()
      {
	if(not isRef)
	  memoryManager<SL>()->release(innerStorage);
      }
#endif
    };
    
    /// Structure to hold statically allocated data
    struct StackStorage
    {
      /// Returns the size
      constexpr Size getSize()
	const
      {
	return
	  StaticSize;
      }
      
      /// Storage
      Fund innerStorage[StaticSize];
      
      /// Return the pointer to inner data
      INLINE_FUNCTION CUDA_HOST_DEVICE
      const Fund* getDataPtr()
	const
      {
	return
	  innerStorage;
      }
      
      PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(getDataPtr,CUDA_HOST_DEVICE);
      
      /// Constructor: since the data is statically allocated, we need to do nothing
      CUDA_HOST_DEVICE
      StackStorage(const Size& size=0)
      {
      }
      
      /// Copy constructor
      CUDA_HOST_DEVICE
      explicit StackStorage(const StackStorage& oth)
      {
	memcpy(this->data,oth.data,StaticSize);
      }
      
      // /// Move constructor is deleted
      // CUDA_HOST_DEVICE
      // StackStorage(StackStorage&&) =delete;
    };
    
    /// Threshold beyond which allocate dynamically in any case
    static constexpr
    Size MAX_STACK_SIZE=2304;
    
    /// Decide whether to allocate on the stack or dynamically
    static constexpr
    bool stackAllocated=
      (StaticSize!=DYNAMIC) and
      (IsStackable==Stackable::MIGHT_GO_ON_STACK) and
      (StaticSize*sizeof(Fund)<=MAX_STACK_SIZE) and
      ((CompilingForDevice==true  and SL==StorLoc::ON_GPU) or
       (CompilingForDevice==false and SL==StorLoc::ON_CPU));
    
    /// Actual storage class
    using ActualStorage=
      std::conditional_t<stackAllocated,StackStorage,DynamicStorage>;
    
    /// Storage of data
    ActualStorage data;
    
    /// Returns the pointer to data
    CUDA_HOST_DEVICE
    Fund* getDataPtr() const
    {
      return data.getDataPtr();
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(getDataPtr,CUDA_HOST_DEVICE);
    
    /// Returns the size
    constexpr Size getSize()
      const
    {
      return
	data.getSize();
    }
    
    /// Allocate the given size
    void allocate(const Size& extSize)
    {
      data.allocate(extSize);
    }
    
    /// Construct taking the size to allocate
    TensorStorage(const Size& size) ///< Size to allocate
      : data(size)
    {
    }
    
    /// Creates starting from a reference
    CUDA_HOST_DEVICE
    TensorStorage(Fund* oth,
		const Size& size) :
      data(oth,size)
    {
      static_assert(stackAllocated==false,"Only dynamic allocation is possible when creating a reference");
    }
    
    /// Default constructor
    CUDA_HOST_DEVICE
    TensorStorage()
    {
    }
    
    /// Copy constructor
    CUDA_HOST_DEVICE
    TensorStorage(const TensorStorage& oth) :
      data(oth.data)
    {
    }
    
    /// Move constructor
    CUDA_HOST_DEVICE
    TensorStorage(TensorStorage&& oth) :
      data(std::move(oth.data))
    {
    }
    
    /// Copy assignemnt
    CUDA_HOST_DEVICE
    TensorStorage& operator=(const TensorStorage& oth)
    {
      data=oth.data;
      
      return *this;
    }
    
    /// Move assignemnt
    TensorStorage& operator=(TensorStorage&& oth)
    {
      data=std::move(oth.data);
      
      return *this;
    }
    
    /// Access to a sepcific value via subscribe operator
    template <typename T>                  // Subscribed component type
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    const Fund& operator[](const T& t)   ///< Subscribed component
      const
    {
      return data.innerStorage[t];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(operator[],CUDA_HOST_DEVICE);
  };
}

#endif
