#ifndef _LOCAL_REDUCTION_V1_HPP
#define _LOCAL_REDUCTION_V1_HPP

#include "base/memory_manager.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  /// Original destructive pairwise reduction, retained for numerical and performance comparison
  template <typename T,
	    typename B,
	    typename Op>
  void locReduceV1(T* loc_res,
		 B& buf,
		 int64_t n,
		 const int nslices,
		 Op&& op)
  {
    if(n%nslices)
      CRASH("number of elements %ld not divisible by number of slices %d",n,nslices);
    
    const int64_t nori_per_slice=n/nslices;
    int64_t nper_slice=n/nslices;
    VERBOSITY_LV3_MASTER_PRINTF("n: %ld, nslices: %d, nori_per_slice: %ld nper_slice: %ld\n",n,nslices,nori_per_slice,nper_slice);
    
    // const double init_time=take_time();
    while(nper_slice>1)
      {
	const int64_t stride=(nper_slice+1)/2;
	const int64_t nreductions_per_slice=nper_slice/2;
	const int64_t nreductions=nreductions_per_slice*nslices;
	VERBOSITY_LV3_MASTER_PRINTF("nper_slice: %ld, stride: %ld, nreductions_per_slice: %ld, nreductions: %ld \n",nper_slice,stride,nreductions_per_slice,nreductions);
	
	PAR(0,
	    nreductions,
	    CAPTURE(nslices,
		    stride,
		    nori_per_slice,
		    op,
		    TO_WRITE(buf)),
	    ireduction,
	  {
	    const int64_t islice=ireduction%nslices;
	    const int64_t ireduction_in_slice=ireduction/nslices;
	    const int64_t first=ireduction_in_slice+nori_per_slice*islice;
	    const int64_t second=first+stride;
	    
	    op(buf[first],buf[second]);
	  });
	  
	  nper_slice=stride;
      }
    
    // MASTER_PRINTF("reduction ended, took %lg s\n",take_time()-init_time);
    
    for(int islice=0;islice<nslices;islice++)
      {
	if constexpr(std::is_standard_layout_v<T> and std::is_trivial_v<T> and not std::is_array_v<T>)
	  {
#ifdef USE_CUDA
	    cudaMemcpy(&loc_res[islice],buf.template getPtrTo<defaultMemorySpace>(islice*nori_per_slice,0),sizeof(T),cudaMemcpyDeviceToHost);
#else
	    loc_res[islice]=buf[islice*nori_per_slice];
#endif
	  }
	else
	  for(size_t idof=0;idof<std::extent_v<T>;idof++)
	    {
#ifdef USE_CUDA
	      cudaMemcpy(&loc_res[islice][idof],buf.template getPtrTo<defaultMemorySpace>(islice*nori_per_slice,idof),sizeof(std::remove_extent_t<T>),cudaMemcpyDeviceToHost);
#else
	      loc_res[islice][idof]=buf(islice*nori_per_slice,idof);
#endif
	    }
      }
  }
  
}

#endif
