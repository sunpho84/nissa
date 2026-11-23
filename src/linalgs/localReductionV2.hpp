#ifndef _LOCAL_REDUCTION_V2_HPP
#define _LOCAL_REDUCTION_V2_HPP

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include "localReductionV1.hpp"

/// V2 reduces each logical slice through independent partials and a final merge
/// Readers separate the algorithm from field layout and from compact scratch storage
/// CPU partials are independent chunks; CUDA threads first accumulate in registers,
/// then merge through block-local shared memory and write one partial per block
/// A second CUDA launch merges those partials, followed by one compact host transfer
/// Every accumulator starts from input, preserving negative maxima without an identity
/// The input is read-only, so tuning can repeat both kernels without backup or refill
/// CUDA uses runtime block dimensions and dynamically sized shared memory by default
/// Unsupported value representations retain the explicitly reported v1 fallback
namespace nissa
{
  namespace localReductionDetails
  {
    /// Copyable reduction value, including C-array components
    template <typename T>
    struct Value
    {
      /// Native value passed to the existing reduction functor
      T data;
    };
    
    /// Copy scalar or array values without flattening across array subobjects
    template <typename T>
    CUDA_HOST_AND_DEVICE inline void assign(T& destination,
                                           const T& source)
    {
      if constexpr(std::is_array_v<T>)
        for(size_t i=0;i<std::extent_v<T>;i++)
          assign(destination[i],source[i]);
      else
        destination=source;
    }
    
    /// Load one value through the field's layout-aware flattened component accessor
    template <typename T,
              typename View>
    CUDA_HOST_AND_DEVICE inline void load(T& destination,
                                         const View& source,
                                         const int64_t site,
                                         int& component)
    {
      if constexpr(std::is_array_v<T>)
        for(size_t i=0;i<std::extent_v<T>;i++)
          load(destination[i],source,site,component);
      else
        destination=source(site,component++);
    }
    
    /// Read original field values without exposing or assuming their physical layout
    template <typename T,
              typename View>
    struct FieldReader
    {
      /// Lightweight readable field reference
      View view;
      
      /// Materialize the native scalar or array for the reduction functor
      CUDA_HOST_AND_DEVICE Value<T> operator()(const int64_t site) const
      {
        Value<T> value;
        int component=0;
        load(value.data,view,site,component);
        
        return value;
      }
    };
    
    /// Read compact partial values in the second reduction stage
    template <typename T>
    struct PartialReader
    {
      /// Contiguous partials, grouped by slice
      const Value<T>* data;
      
      /// Return one already materialized partial
      CUDA_HOST_AND_DEVICE Value<T> operator()(const int64_t site) const
      {
        return data[site];
      }
    };
    
    /// CPU partitioning reference: blockSize*8 gives a target of 2048 values per partial
    /// This is not a CPU thread count and does not select the CUDA block size
    /// CUDA uses the locally declared blockSize=blockDim.x inside reduceBlocks
    inline constexpr int blockSize=256;
    
    /// Maximum partial values per slice, used by both the CPU and CUDA paths
    /// This caps temporary storage and final-stage work, not threads or input length
    /// Above this cap, CPU chunks grow and CUDA threads visit more strided values
    inline constexpr int maxBlocksPerSlice=1024;
    
    /// Preserve the original path for value representations unsupported by the fast implementation
    template <typename T,
              typename B,
              typename Op>
    void fallback(T* result,
                  B& buffer,
                  const int64_t n,
                  const int nslices,
                  const Op& op)
    {
      static const bool reported=[]()
      {
        MASTER_PRINTF("Local reduction v2: using v1 for this value type (size %zu)\n",sizeof(T));
        
        return true;
      }();
      (void)reported;
      locReduceV1(result,buffer,n,nslices,op);
    }
    
#ifdef USE_CUDA
    /// Reduce strided thread accumulators within a block, without assuming a neutral value
    template <typename T,
              typename Reader,
              typename Op>
    __global__ void reduceBlocks(const Reader source,
                                 const int64_t perSlice,
                                 const int blocksPerSlice,
                                 Value<T>* destination,
                                 const Op op)
    {
      // The third launch parameter supplies blockDim.x*sizeof(Value<T>) bytes per block
      extern __shared__ __align__(16) unsigned char storage[];
      auto* partial=reinterpret_cast<Value<T>*>(storage);
      const int blockSize=blockDim.x;
      const int thread=threadIdx.x;
      
      // Blocks are grouped by slice; each writes a separate scratch entry
      const int64_t slice=blockIdx.x/blocksPerSlice;
      const int block=blockIdx.x%blocksPerSlice;
      const int64_t first=static_cast<int64_t>(block)*blockSize+thread;
      
      // Only the first active lanes own values, especially in short final slices
      const int64_t remaining=perSlice-static_cast<int64_t>(block)*blockSize;
      const int active=remaining<blockSize?static_cast<int>(remaining):blockSize;
      if(first<perSlice)
        {
          // Seed from real input, then cover the rest with disjoint grid-stride sequences
          auto value=source(slice*perSlice+first);
          for(int64_t site=first+static_cast<int64_t>(blocksPerSlice)*blockSize;
              site<perSlice;site+=static_cast<int64_t>(blocksPerSlice)*blockSize)
            {
              const auto next=source(slice*perSlice+site);
              op(value.data,next.data);
            }
          partial[thread]=value;
        }
      
      // Inactive lanes must join every barrier, but their uninitialized slots are never read
      __syncthreads();
      for(int stride=blockSize/2;stride>0;stride/=2)
        {
          if(thread<stride and thread+stride<active)
            op(partial[thread].data,partial[thread+stride].data);
          __syncthreads();
        }
      
      // One value per block is read by the second launch, which applies the same algorithm
      if(thread==0)
        destination[blockIdx.x]=partial[0];
    }
    
    /// Aim for eight input values per CUDA thread, with at most maxBlocksPerSlice partials
    /// The grid-stride loop still consumes every input when this heuristic caps the grid
    inline int blocksPerSlice(const int64_t perSlice,
                              const int threads)
    {
      return static_cast<int>(std::min<int64_t>((perSlice-1)/(threads*8)+1,maxBlocksPerSlice));
    }
    
    /// Register each type/layout/operator and slice count separately in the shared tuning archive
    template <typename Reader,
              typename Op>
    int tuningKernelId(const int slices)
    {
      static std::map<int,int> ids;
      auto [entry,inserted]=ids.try_emplace(slices);
      if(inserted)
        {
          const std::string name=
            std::string("localReductionV2:dynamic-v1:")+typeid(Reader).name()+":"+
            typeid(Op).name()+":slices="+std::to_string(slices);
          size_t id=0;
          while(id<kernelInfoLaunchParsStats.size() and kernelInfoLaunchParsStats[id].info.name!=name)
            id++;
          
          if(id==kernelInfoLaunchParsStats.size())
            kernelInfoLaunchParsStats.emplace_back(name);
          
          entry->second=static_cast<int>(id);
        }
      
      return entry->second;
    }
    
    /// Tune read-only reductions without executing backup hooks left by another parallel loop
    template <typename Launch>
    int tuneBlockSize(const int id,
                      const int64_t perSlice,
                      const int maxThreads,
                      const Launch& launch)
    {
      auto beginActions=std::move(benchmarkBeginActions);
      auto endActions=std::move(benchmarkEndActions);
      benchmarkBeginActions.clear();
      benchmarkEndActions.clear();
      const int threads=getOptimalBlockSize(id,perSlice,32,maxThreads,launch);
      benchmarkBeginActions=std::move(beginActions);
      benchmarkEndActions=std::move(endActions);
      if(threads<32 or threads>maxThreads or (threads&(threads-1))!=0)
        CRASH("Invalid cached local reduction block size %d (maximum %d)",threads,maxThreads);
      
      return threads;
    }
#endif
  }
  
  /// Blocked arithmetic reduction; all blocks are seeded from actual input values
  template <typename T,
            typename B,
            typename Op>
  void locReduceV2Arithmetic(T* result,
			     B& buffer,
			     const int64_t n,
			     const int nslices,
			     const Op& op)
  {
    using namespace localReductionDetails;
    const int64_t perSlice=n/nslices;
#ifdef USE_CUDA
    // Limit candidate sizes using the device and both kernels' resource limits
    if constexpr(std::is_arithmetic_v<std::remove_all_extents_t<T>> and sizeof(T)<=128)
      {
        const auto source=buffer.getReadable();
        const FieldReader<T,std::remove_cvref_t<decltype(source)>> reader{source};
        using Reader=std::remove_cvref_t<decltype(reader)>;
        static const int maxThreads=[]()
        {
          cudaFuncAttributes first{},last{};
          DECRYPT_CUDA_ERROR(cudaFuncGetAttributes(&first,reduceBlocks<T,Reader,Op>),"Reduction partial attributes");
          DECRYPT_CUDA_ERROR(cudaFuncGetAttributes(&last,reduceBlocks<T,PartialReader<T>,Op>),"Reduction final attributes");
          int device=0;
          cudaDeviceProp properties{};
          DECRYPT_CUDA_ERROR(cudaGetDevice(&device),"Reduction device");
          DECRYPT_CUDA_ERROR(cudaGetDeviceProperties(&properties,device),"Reduction device properties");
          const size_t sharedAvailable=properties.sharedMemPerBlock-std::max(first.sharedSizeBytes,last.sharedSizeBytes);
          const int limit=std::min({first.maxThreadsPerBlock,last.maxThreadsPerBlock,
                                    static_cast<int>(sharedAvailable/sizeof(Value<T>))});
          if(limit<64)
            CRASH("Insufficient resources for local reduction tuning");
          
          return limit;
        }();
        
        // The smallest candidate (32 threads) needs the most partials
        // Reserve that capacity once; each candidate overwrites its own partials and final output
        const int64_t nPartials=static_cast<int64_t>(nslices)*blocksPerSlice(perSlice,32);
        if(nPartials>std::numeric_limits<int>::max() or
           nPartials+nslices>std::numeric_limits<Size>::max()/sizeof(Value<T>))
          CRASH("Local reduction temporary or launch dimensions overflow");
        
        auto* manager=memoryManager<MemorySpace::GPU>();
        auto* partials=manager->template provide<Value<T>>(nPartials+nslices);
        auto* compact=partials+nPartials;
        const auto launch=[&](const int threads)
        {
          const int blocks=blocksPerSlice(perSlice,threads);
          
          // Dynamic shared memory scales with this candidate, without a separate kernel specialization
          const size_t sharedBytes=threads*sizeof(Value<T>);
          reduceBlocks<T><<<static_cast<unsigned>(static_cast<int64_t>(nslices)*blocks),threads,sharedBytes>>>
            (reader,perSlice,blocks,partials,op);
          DECRYPT_CUDA_ERROR(cudaGetLastError(),"Local reduction v2 partial kernel");
          
          // The compact reader merges only this candidate's initialized partial entries
          reduceBlocks<T><<<nslices,threads,sharedBytes>>>(PartialReader<T>{partials},blocks,1,compact,op);
          DECRYPT_CUDA_ERROR(cudaGetLastError(),"Local reduction v2 final kernel");
        };
        const auto measuredLaunch=[&](const int threads)
        {
          launch(threads);
          
          // Synchronize once after both launches so tuning measures completed device work
          DECRYPT_CUDA_ERROR(cudaDeviceSynchronize(),"Local reduction v2 tuning");
        };
        const int threads=tuneBlockSize(tuningKernelId<Reader,Op>(nslices),perSlice,maxThreads,measuredLaunch);
        launch(threads);
        std::vector<Value<T>> host(nslices);
        DECRYPT_CUDA_ERROR(cudaMemcpy(host.data(),compact,nslices*sizeof(Value<T>),cudaMemcpyDeviceToHost),
                           "Local reduction v2 result transfer");
        for(int slice=0;slice<nslices;slice++)
          assign(result[slice],host[slice].data);
        manager->release(partials);
      }
    else
      fallback(result,buffer,n,nslices,op);
#else
    const int blocks=static_cast<int>(std::min<int64_t>((perSlice-1)/(blockSize*8)+1,maxBlocksPerSlice));
    const int64_t nPartials=static_cast<int64_t>(nslices)*blocks;
    const auto source=buffer.getReadable();
    const FieldReader<T,std::remove_cvref_t<decltype(source)>> reader{source};
    std::vector<Value<T>> partials(nPartials);
#ifdef HAVE_OPENMP
#pragma omp parallel for if(nPartials>1)
#endif
    for(int64_t partial=0;partial<nPartials;partial++)
      {
        const int64_t slice=partial/blocks;
        const int block=partial%blocks;
        const int64_t chunk=perSlice/blocks;
        const int64_t remainder=perSlice%blocks;
        const int64_t begin=block*chunk+std::min<int64_t>(block,remainder);
        const int64_t end=begin+chunk+(block<remainder);
        auto value=reader(slice*perSlice+begin);
        for(int64_t site=begin+1;site<end;site++)
          {
            const auto next=reader(slice*perSlice+site);
            op(value.data,next.data);
          }
        partials[partial]=value;
      }
    for(int slice=0;slice<nslices;slice++)
      {
        auto value=partials[static_cast<int64_t>(slice)*blocks];
        for(int block=1;block<blocks;block++)
          op(value.data,partials[static_cast<int64_t>(slice)*blocks+block].data);
        assign(result[slice],value.data);
      }
#endif
  }
  
  /// Select the fast arithmetic path without imposing new type requirements on v1 callers
  template <typename T,
            typename B,
            typename Op>
  void locReduceV2(T* result,
                   B& buffer,
                   const int64_t n,
                   const int nslices,
                   const Op& op)
  {
    if constexpr(std::is_arithmetic_v<std::remove_all_extents_t<T>>)
      locReduceV2Arithmetic(result,buffer,n,nslices,op);
    else
      localReductionDetails::fallback(result,buffer,n,nslices,op);
  }

}

#endif
