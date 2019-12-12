#ifndef _MEMORY_MANAGER_HPP
#define _MEMORY_MANAGER_HPP

#include <map>
#include <vector>
#include <cstdint>

#ifndef EXTERN_MEMORY_MANAGER
 #define EXTERN_MEMORY_MANAGER extern
#endif

#include <base/debug.hpp>
#include <routines/ios.hpp>
#include <new_types/value_with_extreme.hpp>

namespace nissa
{
  /// Type used for size
  using Size=int64_t;
  
  /// Minimal alignment
#define DEFAULT_ALIGNMENT 16
  
  /// Memory manager
  class MemoryManager
  {
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
    
    /// Number of unaligned allocation performed
    Size nUnalignedAlloc{0};
    
    /// Number of aligned allocation performed
    Size nAlignedAlloc{0};
    
    /// Number of cached memory reused
    Size nCachedReused{0};
    
    /// Get aligned memory
    ///
    /// Call the system routine which allocate memory
    void* allocateRawAligned(const Size size,         ///< Amount of memory to allocate
			     const Size alignment)    ///< Required alignment
    {
      // runLog()<<"Raw allocating "<<size;
      
      /// Result
      void* ptr=nullptr;
      
      /// Returned condition
      int rc=posix_memalign(&ptr,alignment,size);
      
      if(rc) crash("Failed to allocate %ld with alignement %ld",size,alignment);
      
      nAlignedAlloc++;
      
      return ptr;
    }
    
    /// Add to the list of used memory
    void pushToUsed(void* ptr,
		    const Size size)
    {
      used[ptr]=size;
      
      usedSize+=size;
      
      // runLog()<<"Pushing to used "<<ptr<<" "<<size<<", number of used:"<<used.size();
    }
    
    /// Removes a pointer from the used list, without actually freeing associated memory
    ///
    /// Returns the size of the memory pointed
    Size popFromUsed(void* ptr) ///< Pointer to the memory to move to cache
    {
      // runLog()<<"Popping from used "<<ptr;
      
      /// Iterator to search result
      auto el=used.find(ptr);
      
      if(el==used.end())
	crash("Unable to find dinamically allocated memory %p",ptr);
      
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
      
      // runLog()<<"Pushing to cache "<<size<<" "<<ptr<<", cache size: "<<cached.size();
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
      // runLog()<<"Try to popping from cache "<<size;
      
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
      // runLog()<<"Moving to cache "<<ptr;
      
      /// Size of pointed memory
      const Size size=popFromUsed(ptr);
      
      pushToCache(ptr,size);
    }
    
  public:
    
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
    
    /// Allocate or get from cache after computing the proper size
    template <class T>
    T* provideAligned(const Size nel,
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
	ptr=allocateRawAligned(size,alignment);
      else
	nCachedReused++;
      
      pushToUsed(ptr,size);
      
      return static_cast<T*>(ptr);
    }
    
    /// Decleare unused the memory
    template <class T>
    void release(T* ptr) ///< Pointer getting freed
    {
      if(useCache)
	moveToCache(static_cast<void*>(ptr));
      else
	{
	  popFromUsed(ptr);
	  free(ptr);
	}
    }
    
    /// Release all used memory
    void releaseAllUsedMemory()
    {
      /// Iterator on elements to release
      auto el=
	used.begin();
      
      while(el!=used.end())
	{
	  // runLog()<<"Releasing "<<el.first<<" size "<<el.second;
	  
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
	      // runLog()<<"Removing from cache size "<<el.first;
	      
	      /// Memory to free
	      void* ptr=popFromCache(size,DEFAULT_ALIGNMENT);
	      
	      free(ptr);
	    }
	}
    }
    
    /// Print to a stream
    void printStatistics()
    {
      master_printf("Maximal memory used: %ld bytes, currently used: %ld bytes, number of allocations: %ld unaligned, %ld aligned\n",
		    usedSize.extreme(),(Size)usedSize,nUnalignedAlloc,nAlignedAlloc);
      master_printf("Maximal memory cached: %ld bytes, currently used: %ld bytes, number of reused: %ld\n",
		    cachedSize.extreme(),(Size)cachedSize,nCachedReused);
    }
    
    /// Create the memory manager
    MemoryManager() :
      usedSize(0),
      cachedSize(0)
    {
      master_printf("Starting the memory manager\n");
    }
    
    /// Destruct the memory manager
    ~MemoryManager()
    {
      master_printf("Stopping the memory manager\n");
      
      printStatistics();
      
      releaseAllUsedMemory();
      
      clearCache();
    }
  };
  
  EXTERN_MEMORY_MANAGER MemoryManager *memory_manager;
}

#undef EXTERN_MEMORY_MANAGER

#endif
