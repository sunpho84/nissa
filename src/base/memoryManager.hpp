#ifndef _MEMORY_MANAGER_HPP
#define _MEMORY_MANAGER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <map>
#include <vector>
#include <cstdint>

#ifdef ENABLE_DEVICE_CODE
# include <cuda_runtime.h>
#endif

#ifdef ENABLE_DEVICE_CODE
# include <base/cuda.hpp>
#endif
#include <base/memoryType.hpp>
#include <expr/execSpace.hpp>
#include <newTypes/valueWithExtreme.hpp>
#include <metaprogramming/constnessChanger.hpp>
#include <metaprogramming/crtp.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  /// Class needed to dispatch not allocation
  struct DoNotAllocate
  {
  };
  
  /// Handle to avoid allocation
  constexpr DoNotAllocate doNotAllocate;
  
  /////////////////////////////////////////////////////////////////
  
  /// Type used for size
  using Size=int64_t;
  
  /// Minimal alignment
#define DEFAULT_ALIGNMENT 16
  
  /// Memory manager, base type
  struct MemoryManager
  {
    /// Number of allocation performed
    Size nAlloc{0};
    
    /// List of dynamical allocated memory
    std::map<void*,Size> used;
    
    /// List of cached allocated memory
    std::map<Size,std::vector<void*>> cached;
    
    /// Size of used memory
    ValWithMax<Size> usedSize;
    
    /// Size of cached memory
    ValWithMax<Size> cachedSize;
    
    /// Use or not cache
    bool useCache{true};
    
    /// Number of cached memory reused
    Size nCachedReused{0};
    
    /// Add to the list of used memory
    void pushToUsed(void* ptr,
		    const Size size)
    {
      used[ptr]=size;
      
      usedSize+=size;
      
      VERBOSITY_LV3_MASTER_PRINTF("Pushing to used %p %zu, used: %zu\n",ptr,size,used.size());
    }
    
    /// Removes a pointer from the used list, without actually freeing associated memory
    ///
    /// Returns the size of the memory pointed
    Size popFromUsed(void* ptr) ///< Pointer to the memory to move to cache
    {
      VERBOSITY_LV3_MASTER_PRINTF("Popping from used %p\n",ptr);
      
      /// Iterator to search result
      auto el=used.find(ptr);
      
      if(el==used.end())
	CRASH("Unable to find dinamically allocated memory %p",ptr);
      
      /// Size of memory
      const Size size=el->second;
      
      usedSize-=size;
      
      used.erase(el);
      
      return size;
    }
    
    /// Adds a memory to cache
    void pushToCache(void* ptr,          ///< Memory to cache
		     const Size size)    ///< Memory size
    {
      cached[size].push_back(ptr);
      
      cachedSize+=size;
      
      VERBOSITY_LV3_MASTER_PRINTF("Pushing to cache %zu %p, cache size: %zu\n",size,ptr,cached.size());
    }
    
    /// Check if a pointer is suitably aligned
    static bool isAligned(const void* ptr,
			  const Size alignment)
    {
      return reinterpret_cast<uintptr_t>(ptr)%alignment==0;
    }
    
    /// Pop from the cache, returning to use
    void* popFromCache(const Size& size,
		       const Size& alignment)
    {
      VERBOSITY_LV3_MASTER_PRINTF("Try to popping from cache %zu\n",size);
      
      /// List of memory with searched size
      auto cachedIt=cached.find(size);
      
      if(cachedIt==cached.end())
	return nullptr;
      else
	{
	  /// Vector of pointers
	  auto& list=cachedIt->second;
	  
	  /// Get latest cached memory with appropriate alignment
	  auto it=list.end()-1;
	  
	  while(it!=list.begin()-1 and not isAligned(*it,alignment))
	    it--;
	  
	  if(it==list.begin()-1)
	    return nullptr;
	  else
	    {
	      /// Returned pointer, copied here before erasing
	      void* ptr=*it;
	      
	      list.erase(it);
	      
	      cachedSize-=size;
	      
	      if(list.size()==0)
		cached.erase(cachedIt);
	      
	      return ptr;
	    }
	}
    }
    
    /// Move the allocated memory to cache
    void moveToCache(void* ptr) ///< Pointer to the memory to move to cache
    {
      VERBOSITY_LV3_MASTER_PRINTF("Moving to cache %p\n",ptr);
      
      /// Size of pointed memory
      const Size size=popFromUsed(ptr);
      
      pushToCache(ptr,size);
    }
    
    /// Enable cache usage
    void enableCache()
    {
      useCache=true;
    }
    
    /// Disable cache usage
    void disableCache()
    {
      useCache=false;
      
      clearCache();
    }
    
    virtual void* allocateRaw(const Size& size,            ///< Amount of memory to allocate
			      const Size& alignment) =0;   ///< Required alignment
    
    /// Allocate or get from cache after computing the proper size
    template <class T>
    T* provide(const Size nel,
	       const Size alignment=DEFAULT_ALIGNMENT)
    {
      /// Total size to allocate
      const Size size=sizeof(T)*nel;
      
      /// Allocated memory
      void* ptr;
      
      // Search in the cache
      ptr=popFromCache(size,alignment);
      
      // If not found in the cache, allocate new memory
      if(ptr==nullptr)
	ptr=allocateRaw(size,alignment);
      else
	nCachedReused++;
      
      pushToUsed(ptr,size);
      
      return static_cast<T*>(ptr);
    }
    
    /// Properly free
    virtual void deAllocateRaw(void* ptr) =0;
    
    /// Declare unused the memory and possibly free it
    template <typename T>
    void release(T* &ptr,                      ///< Pointer getting freed
		 const bool& forceFree=false)  ///< Force to free and not to move to cache
    {
      if(useCache and not forceFree)
	moveToCache(ptr);
      else
	{
	  popFromUsed(ptr);
	  deAllocateRaw(ptr);
	}
      
      ptr=nullptr;
    }
    
    /// Release all used memory
    void releaseAllUsedMemory()
    {
      VERBOSITY_LV2_MASTER_PRINTF("Releasing all used memory\n");
      
      /// Iterator on elements to release
      auto el=
	used.begin();
      
      while(el!=used.end())
	{
	  VERBOSITY_LV3_MASTER_PRINTF("Releasing %p size %zu\n",el->first,el->second);
	  
	  /// Pointer to memory to release
	  void* ptr=el->first;
	  
	  // Increment iterator before releasing
	  el++;
	  
	  release(ptr);
	}
    }
    
    /// Release all memory from cache
    void clearCache()
    {
      VERBOSITY_LV2_MASTER_PRINTF("Clearing cache\n");
      
      /// Iterator to elements of the cached memory list
      auto el=cached.begin();
      
      while(el!=cached.end())
	{
	  /// Number of elements to free
	  const Size n=el->second.size();
	  
	  /// Size to be removed
	  const Size size=el->first;
	  
	  // Increment before erasing
	  el++;
	  
	  for(Size i=0;i<n;i++)
	    {
	      VERBOSITY_LV3_MASTER_PRINTF("Removing from cache size %zu ",el->first);
	      
	      /// Memory to free
	      void* ptr=popFromCache(size,DEFAULT_ALIGNMENT);
	      
	      VERBOSITY_LV3_MASTER_PRINTF("ptr: %p\n",ptr);
	      deAllocateRaw(ptr);
	    }
	}
    }
    
    /// Print to a stream
    void printStatistics()
    {
      masterPrintf("Maximal memory used: %ld bytes, currently used: %ld bytes, number of allocation: %ld\n",
		   usedSize.extreme(),(Size)usedSize,nAlloc);
      masterPrintf("Maximal memory cached: %ld bytes, currently used: %ld bytes, number of reused: %ld\n",
		    cachedSize.extreme(),(Size)cachedSize,nCachedReused);
    }
    
    /// Create the memory manager
    MemoryManager() :
      usedSize(0),
      cachedSize(0)
    {
      masterPrintf("Starting the memory manager\n");
    }
    
    /// Base destructor
    virtual ~MemoryManager() =default;
    
    /// Destruct the memory manager
    void destroy()
    {
      masterPrintf("Stopping the memory manager\n");
      
      printStatistics();
      
      disableCache();
      
      releaseAllUsedMemory();
    }
  };
  
  /// Manager of CPU memory
  struct CPUMemoryManager :
    MemoryManager
  {
    /// Get memory
    ///
    /// Call the system routine which allocate memory
    void* allocateRaw(const Size& size,        ///< Amount of memory to allocate
		      const Size& alignment)   ///< Required alignment
    {
      /// Result
      void* ptr=nullptr;
      
      /// Returned condition
      VERBOSITY_LV3_MASTER_PRINTF("Allocating size %zu on CPU, ",size);
      int rc=posix_memalign(&ptr,alignment,size);
      if(rc)
	CRASH("Failed to allocate %ld CPU memory with alignement %ld",size,alignment);
      VERBOSITY_LV3_MASTER_PRINTF("ptr: %p\n",ptr);
      
      nAlloc++;
      
      return ptr;
    }
    
    /// Properly free
    void deAllocateRaw(void* ptr)
    {
      VERBOSITY_LV3_MASTER_PRINTF("Freeing from CPU memory %p\n",ptr);
      free(ptr);
    }
    
    /// Destruct calling common procedure
    ~CPUMemoryManager()
    {
      destroy();
    }
  };
  
  inline MemoryManager* cpuMemoryManager;
  
#ifdef ENABLE_DEVICE_CODE
  
  /// Manager of GPU memory
  struct GPUMemoryManager :
    MemoryManager
  {
    /// Get memory on GPU
    ///
    /// Call the system routine which allocate memory
    void* allocateRaw(const Size& size,        ///< Amount of memory to allocate
		      const Size& alignment)   ///< Required alignment
    {
      /// Result
      void* ptr=nullptr;
      
      VERBOSITY_LV3_MASTER_PRINTF("Allocating size %zu on GPU, ",size);
      decryptCudaError(cudaMalloc(&ptr,size),"Allocating on Gpu");
      VERBOSITY_LV3_MASTER_PRINTF("ptr: %p\n",ptr);
      
      nAlloc++;
      
      return ptr;
    }
    
    /// Properly free
    void deAllocateRaw(void* ptr)
    {
      masterPrintf("Freeing from GPU memory %p\n",ptr);
      decryptCudaError(cudaFree(ptr),"Freeing from GPU");
    }
    
    /// Destruct calling common procedure
    ~GPUMemoryManager()
    {
      destroy();
    }
  };
  
  inline MemoryManager* gpuMemoryManager;
  
#endif
  
  template <MemoryType MT>
  inline MemoryManager* memoryManager()
  {
    switch(MT)
      {
      case MemoryType::CPU:
	return cpuMemoryManager;
	break;
#ifdef ENABLE_DEVICE_CODE
      case MemoryType::GPU:
	return gpuMemoryManager;
	break;
#endif
    }
  }
  
#ifdef ENABLE_DEVICE_CODE
  
  template <MemoryType FROM,
	    MemoryType TO>
  cudaMemcpyKind cudaMemcpyKindForTransferFromTo;
  
#define PROVIDE_CUDA_MEMCPY_KIND_FOR_TRANSFER(MFROM,MTO,CFROM,CTO)	\
  template <>								\
  constexpr								\
  cudaMemcpyKind cudaMemcpyKindForTransferFromTo<MemoryType::MFROM,MemoryType::MTO> =cudaMemcpy ## CFROM ## To ## CTO
  
  PROVIDE_CUDA_MEMCPY_KIND_FOR_TRANSFER(CPU,CPU,Host,Host);
  PROVIDE_CUDA_MEMCPY_KIND_FOR_TRANSFER(CPU,GPU,Host,Device);
  PROVIDE_CUDA_MEMCPY_KIND_FOR_TRANSFER(GPU,CPU,Device,Host);
  PROVIDE_CUDA_MEMCPY_KIND_FOR_TRANSFER(GPU,GPU,Device,Device);
  
#undef PROVIDE_CUDA_MEMCPY_KIND_FOR_TRANSFER
  
#endif
  
  template <MemoryType TO,
	    MemoryType FROM>
  void memcpy(void* dst,
	      const void* src,
	      const size_t& count)
  {
#ifdef ENABLE_DEVICE_CODE
    decryptCudaError(cudaMemcpy(dst,src,count,cudaMemcpyKindForTransferFromTo<FROM,TO>),"calling cudaMemcpy");
#else
    ::memcpy(dst,src,count);
#endif
  }
  
  /// Returns the memory manager suitable for the execution space
  template <ExecSpace ES>
  inline MemoryManager* memoryManager()
  {
    return memoryManager<getMemoryType<ES>()>();
  }
  
  /// Initialize the memory managers
  inline void initMemoryManagers()
  {
    //initialize the memory manager
    cpuMemoryManager=new CPUMemoryManager;
#ifdef ENABLE_DEVICE_CODE
    gpuMemoryManager=new GPUMemoryManager;
#endif
  }
}

#endif
