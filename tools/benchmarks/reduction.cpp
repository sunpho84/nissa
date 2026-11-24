#include "nissa.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstring>
#include <memory>
#include <limits>
#include <string_view>
#include <vector>

namespace
{
  using namespace nissa;
  
  /// Non-arithmetic scalar exercising the compatibility fallback
  struct PairValue
  {
    /// First accumulated component
    double first;
    
    /// Second accumulated component
    double second;
    
    /// Support the original scalar sum functor
    CUDA_HOST_AND_DEVICE PairValue& operator+=(const PairValue& other)
    {
      first+=other.first;
      second+=other.second;
      
      return *this;
    }
  };
  
  /// Deterministic mixed-sign sums and strictly negative maxima, independent of the RNG
  CUDA_HOST_AND_DEVICE double sample(const int64_t site,
                                     const int component,
                                     const int sourceRank,
                                     const bool maximum)
  {
    return maximum?-1.0-((site%97)+component+sourceRank)*0.125:
      ((site%101)-50)*(component+1)*0.125+(sourceRank+1)*0.0001;
  }
  
  /// Populate outside reduction timings, since v1 consumes the input
  template <typename T,
            SpaceTimeLayout layout>
  void fill(Field<T,FULL_SPACE,layout>& field,
            const int64_t n,
            const bool maximum)
  {
    const int sourceRank=rank;
    PAR(0,n,CAPTURE(TO_WRITE(field),sourceRank,maximum),site,
    {
      for(int component=0;component<sizeof(T)/sizeof(double);component++)
        field(site,component)=sample(site,component,sourceRank,maximum);
    });
  }
  
  /// Compare against a host reference with a summation-order-aware absolute tolerance
  template <typename T>
  void check(const T* result,
             const int64_t n,
             const int slices,
             const int sourceRank,
             const bool maximum)
  {
    constexpr int components=sizeof(T)/sizeof(double);
    const int64_t perSlice=n/slices;
    for(int slice=0;slice<slices;slice++)
      for(int component=0;component<components;component++)
        {
          long double expected=maximum?-std::numeric_limits<long double>::infinity():0;
          long double magnitude=0;
          for(int64_t site=slice*perSlice;site<(slice+1)*perSlice;site++)
            {
              const double value=sample(site,component,sourceRank,maximum);
              if(maximum)
                expected=std::max(expected,static_cast<long double>(value));
              else
                expected+=value;
              magnitude+=std::abs(value);
            }
          double observed=0;
          if constexpr(std::is_array_v<T>)
            observed=result[slice][component];
          else
            observed=result[slice];
          const long double tolerance=maximum?0:1e-11L+1e-12L*magnitude;
          if(not std::isfinite(observed) or std::abs(observed-expected)>tolerance)
            CRASH("Reduction mismatch rank %d slice %d component %d: %.17g vs %.21Lg (tolerance %.3Lg)",
                  sourceRank,slice,component,observed,expected,tolerance);
        }
  }
  
  /// Verify v2's read-only contract outside all timings
  template <typename T,
            SpaceTimeLayout layout>
  void checkPreserved(const Field<T,FULL_SPACE,layout>& field,
                      const int64_t n,
                      const bool maximum)
  {
    std::vector<double> host(field.nTotalElements());
#ifdef USE_CUDA
    DECRYPT_CUDA_ERROR(cudaMemcpy(host.data(),field.template getPtr<defaultMemorySpace>(),
                                 host.size()*sizeof(double),cudaMemcpyDeviceToHost),
                       "Reduction comparison input transfer");
#else
    std::memcpy(host.data(),field.template getPtr<defaultMemorySpace>(),host.size()*sizeof(double));
#endif
    for(int64_t site=0;site<n;site++)
      for(int component=0;component<sizeof(T)/sizeof(double);component++)
        if(host[field.index(site,component)]!=sample(site,component,rank,maximum))
          CRASH("Local reduction v2 changed its input at site %ld component %d",site,component);
  }
  
#ifdef USE_CUDA
  /// Check every candidate, including partially occupied blocks and multiple reduction stages
  template <typename T,
            SpaceTimeLayout layout>
  void validateBlockSizes()
  {
    using namespace localReductionDetails;
    Field<T,FULL_SPACE,layout> field("reduction launch comparison");
    constexpr int slices=3;
    for(const int64_t perSlice : {1,3,33,257,4097})
      if(perSlice*slices<=locVol)
        {
          const auto run=[&](const auto& op,
                            const bool maximum)
          {
            fill(field,perSlice*slices,maximum);
            const auto source=field.getReadable();
            const FieldReader<T,std::remove_cvref_t<decltype(source)>> reader{source};
            const int64_t capacity=slices*blocksPerSlice(perSlice,32);
            auto* manager=memoryManager<MemorySpace::GPU>();
            auto* partials=manager->template provide<Value<T>>(capacity+slices);
            auto* compact=partials+capacity;
            std::vector<Value<T>> host(slices);
            std::unique_ptr<T[]> result{new T[slices]};
            for(const int threads : {32,64,128,256,512,1024})
              if(threads<=deviceProp.maxThreadsPerBlock and threads*sizeof(Value<T>)<=deviceProp.sharedMemPerBlock)
                {
                  const int blocks=blocksPerSlice(perSlice,threads);
                  const size_t sharedBytes=threads*sizeof(Value<T>);
                  reduceBlocks<T><<<slices*blocks,threads,sharedBytes>>>(reader,perSlice,blocks,partials,op);
                  DECRYPT_CUDA_ERROR(cudaGetLastError(),"Candidate partial reduction");
                  reduceBlocks<T><<<slices,threads,sharedBytes>>>(PartialReader<T>{partials},blocks,1,compact,op);
                  DECRYPT_CUDA_ERROR(cudaGetLastError(),"Candidate final reduction");
                  DECRYPT_CUDA_ERROR(cudaMemcpy(host.data(),compact,slices*sizeof(Value<T>),cudaMemcpyDeviceToHost),
                                     "Candidate reduction result");
                  for(int slice=0;slice<slices;slice++)
                    assign(result[slice],host[slice].data);
                  check(result.get(),perSlice*slices,slices,rank,maximum);
                }
            checkPreserved(field,perSlice*slices,maximum);
            manager->release(partials);
          };
          run(GlbReduceSumFunctor{},false);
          if constexpr(not std::is_array_v<T>)
            run(GlbReduceMaxFunctor{},true);
        }
  }
#endif
  
  /// Exercise slice boundaries, non-power-of-two lengths, both layouts, and negative maxima
  template <typename T,
            SpaceTimeLayout layout>
  void validate()
  {
    Field<T,FULL_SPACE,layout> field("reduction comparison");
    for(const int64_t n : {int64_t{1},int64_t{3},int64_t{9},int64_t{17},int64_t{257},int64_t{771},int64_t{4097},locVol})
      if(n<=locVol)
        for(const int slices : {1,3})
          if(n%slices==0)
            for(const auto version : {LocalReductionVersion::V1,LocalReductionVersion::V2})
              {
                std::unique_ptr<T[]> result{new T[slices]};
                fill(field,n,false);
                locReduce(result.get(),field,n,slices,GlbReduceSumFunctor{},version);
                check(result.get(),n,slices,rank,false);
                if(version==LocalReductionVersion::V2)
                  checkPreserved(field,n,false);
                if constexpr(not std::is_array_v<T>)
                  {
                    fill(field,n,true);
                    locReduce(result.get(),field,n,slices,GlbReduceMaxFunctor{},version);
                    check(result.get(),n,slices,rank,true);
                    if(version==LocalReductionVersion::V2)
                      checkPreserved(field,n,true);
                  }
              }
  }
  
  /// Verify the unchanged MPI stage and nonzero offsets using the selected environment version
  void validateGlobal()
  {
    int ranks=0;
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    const int slices=locSize[0];
    if(ranks>std::numeric_limits<int>::max()/slices)
      CRASH("Reduction comparison global slice count overflows int");
    Field<complex,FULL_SPACE> field("global reduction comparison");
    fill(field,locVol,false);
    std::unique_ptr<complex[]> result{new complex[slices*ranks]};
    glb_reduce(result.get(),field,locVol,slices*ranks,slices,rank*slices);
    for(int sourceRank=0;sourceRank<ranks;sourceRank++)
      check(result.get()+sourceRank*slices,locVol,slices,sourceRank,false);
  }
  
  /// Unsupported representations retain v1 behavior, including writable-field invalidation
  void validateFallback()
  {
    Field<PairValue,FULL_SPACE> field("reduction fallback");
    PAR(0,3,CAPTURE(TO_WRITE(field)),site,
    {
      field[site]=PairValue{1,-2};
    });
    field._haloIsValid=true;
    PairValue result;
    locReduce(&result,field,3,1,GlbReduceSumFunctor{},LocalReductionVersion::V2);
    if(result.first!=3 or result.second!=-6 or field.haloIsValid())
      CRASH("Local reduction compatibility fallback failed");
  }
  
  /// Print maximum-rank local time, excluding input refill, reference checks, and MPI reduction
  template <typename T>
  void benchmark(const int iterations,
                 const LocalReductionVersion version)
  {
    Field<T,FULL_SPACE> field("reduction benchmark");
    const int slices=locSize[0];
    std::unique_ptr<T[]> result{new T[slices]};
    fill(field,locVol,false);
    MPI_Barrier(MPI_COMM_WORLD);
    const double warmupBegin=take_time();
    locReduce(result.get(),field,locVol,slices,GlbReduceSumFunctor{},version);
    const double warmup=take_time()-warmupBegin;
    check(result.get(),locVol,slices,rank,false);
    double total=0;
    for(int iteration=0;iteration<iterations;iteration++)
      {
        fill(field,locVol,false);
        MPI_Barrier(MPI_COMM_WORLD);
        const double begin=take_time();
        locReduce(result.get(),field,locVol,slices,GlbReduceSumFunctor{},version);
        total+=take_time()-begin;
      }
    double maximumTotal=0,maximumWarmup=0;
    MPI_Allreduce(&total,&maximumTotal,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&warmup,&maximumWarmup,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    check(result.get(),locVol,slices,rank,false);
    MASTER_PRINTF("REDUCTION v%d components=%zu local_sites=%ld slices=%d iterations=%d warmup_seconds=%.9g mean_seconds=%.9g\n",
                  static_cast<int>(version),sizeof(T)/sizeof(double),locVol,slices,iterations,
                  maximumWarmup,maximumTotal/iterations);
  }
  
  /// Parse benchmark sizes without accepting partial or out-of-range arguments
  int positiveArgument(const char* text)
  {
    const std::string_view value(text);
    int result=0;
    const auto parsed=std::from_chars(value.data(),value.data()+value.size(),result);
    if(parsed.ec!=std::errc{} or parsed.ptr!=value.data()+value.size() or result<=0)
      CRASH("Expected a positive integer, received '%s'",text);
    
    return result;
  }
}

/// Validate both implementations, then benchmark them with identical refilled input
int main(int narg,
         char** arg)
{
  using namespace nissa;
  initNissa(narg,arg);
  if(narg!=1 and narg!=4)
    CRASH("Usage: reduction [L T iterations]");
  
  const int L=narg==4?positiveArgument(arg[1]):8;
  const int T=narg==4?positiveArgument(arg[2]):16;
  const int iterations=narg==4?positiveArgument(arg[3]):5;
  localReductionVersion();
  initGrid(T,L);
  MASTER_PRINTF("REDUCTION validation begin (v1 tuning may occur during validation)\n");
#ifdef USE_CUDA
  validateBlockSizes<double,SpaceTimeLayout::CPU>();
  validateBlockSizes<double,SpaceTimeLayout::GPU>();
  validateBlockSizes<complex,SpaceTimeLayout::CPU>();
  validateBlockSizes<complex,SpaceTimeLayout::GPU>();
  MASTER_PRINTF("REDUCTION CUDA block-size candidates passed\n");
#endif
  validate<double,SpaceTimeLayout::CPU>();
  validate<double,SpaceTimeLayout::GPU>();
  validate<complex,SpaceTimeLayout::CPU>();
  validate<complex,SpaceTimeLayout::GPU>();
  validateFallback();
  validateGlobal();
  MPI_Barrier(MPI_COMM_WORLD);
  MASTER_PRINTF("REDUCTION validation passed on all ranks\n");
  for(const auto version : {LocalReductionVersion::V1,LocalReductionVersion::V2})
    {
      benchmark<double>(iterations,version);
      benchmark<complex>(iterations,version);
    }
  closeNissa();
  
  return 0;
}
