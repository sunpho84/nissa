#ifndef _REDUCE_HPP
#define _REDUCE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cstdint>

#ifdef USE_CUDA
# include <thrust/complex.h>
# include <thrust/execution_policy.h>
# include <thrust/reduce.h>
#endif

#include "base/memory_manager.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#ifndef EXTERN_REDUCE
# define EXTERN_REDUCE extern
# define INIT_REDUCE_TO(...)
#else
# define INIT_REDUCE_TO(...) __VA_ARGS__
#endif

namespace nissa
{
//   /// Size of the reducing buffer
//   EXTERN_REDUCE int64_t reducing_buffer_size INIT_REDUCE_TO(=0);
  
//   /// Reference to the reducing buffer
//   EXTERN_REDUCE CUDA_MANAGED void *reducing_buffer INIT_REDUCE_TO(=nullptr);
  
//   /// Get the reduction buffer
//   template <typename T>
//   T* get_reducing_buffer(const int64_t& n)
//   {
//     const int64_t required_size=sizeof(T)*n;
//     if(reducing_buffer_size<required_size)
//       {
// 	if(reducing_buffer!=nullptr)
// 	  nissa_free(reducing_buffer);
	
// 	//reducing_buffer_size=0; check
//       }
    
//     if(reducing_buffer==nullptr)
//       {
// 	if(reducing_buffer_size==0)
// 	  MASTER_PRINTF("Allocating the reduction buffer to %ld bytes\n",required_size);
// 	else
// 	  MASTER_PRINTF("Reallocating the reduction buffer from %ld bytes to %ld bytes\n",reducing_buffer_size,required_size);
	
// 	reducing_buffer=nissa_malloc("reducing_buffer",required_size,char);
// 	reducing_buffer_size=required_size;
//       }
    
//     return (T*)reducing_buffer;
//   }
  
//   void deallocate_reduction_buffer();
  
//   void loc_reduce(int64_t* loc_res,int64_t* buf,int64_t n,int nslices=1);
//   void loc_reduce(double* loc_res,double* buf,int64_t n,int nslices=1);
//   void loc_reduce(complex* loc_res,complex* buf,int64_t n,int nslices=1);
// #ifdef REPRODUCIBLE_RUN
//   void loc_reduce(float_128* loc_res,float_128* buf,int64_t n,int nslices=1);
//   void loc_reduce(complex_128* loc_res,complex_128* buf,int64_t n,int nslices=1);
// #endif
  
  //reduce a vector
  template <typename T>
  inline void non_loc_reduce(T* out_glb,T* in_loc=nullptr,const int n=1,MPI_Op mpi_op=MPI_Op_sum_for_type<T>())
  {
    if(in_loc==nullptr)
      in_loc=(T*)MPI_IN_PLACE;
    
    MPI_Allreduce(in_loc,out_glb,n,MPI_Datatype_of<T>(),mpi_op,MPI_COMM_WORLD);
  }
  
  /////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T,
	    typename B,
	    typename Op>
  void locReduce(T* loc_res,
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
	if constexpr(std::is_pod_v<T> and not std::is_array_v<T>)
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
  
  /////////////////////////////////////////////////////////////////
  
  /// Functor to handle the pairwise sum reduce
  struct GlbReduceSumFunctor
  {
    ///Pairwise reduction
    template <typename F1,
	      typename F2>
    CUDA_DEVICE INLINE_FUNCTION
    void operator()(F1&& res,
		    const F2& acc) const
    {
      if constexpr(not std::is_array_v<std::remove_reference_t<F1>>)
	res+=acc;
      else
	for(size_t i=0;i<std::extent_v<std::remove_reference_t<F1>>;i++)
	  GlbReduceSumFunctor::operator()(res[i],acc[i]);
    }
  };
  
  /// Functor to handle the pairwise max reduce
  struct GlbReduceMaxFunctor
  {
    ///Pairwise reduction
    template <typename F1,
	      typename F2>
    CUDA_DEVICE INLINE_FUNCTION
    void operator()(F1&& res,
		    const F2& acc) const
    {
      if(acc>res)
	res=acc;
    }
  };
  
  /// Reduce a vector over all nodes
  template <typename T,
	    typename B>
  void glb_reduce(T* glb_res,
		  B& buf,
		  const int64_t nloc,
		  const int nslices=1,
		  const int nloc_slices=1,
		  const int loc_offset=0)
  {
    /// to be improved, catching fields
    T loc_res[nslices];
    memset(loc_res,0,sizeof(T)*nslices);
    
    locReduce(loc_res+loc_offset,
	      buf,
	      nloc,
	      nloc_slices,
	      GlbReduceSumFunctor());
    
    non_loc_reduce(glb_res,loc_res,nslices);
  }
  
  /// Reduce a vector over all nodes
  template <typename T,
	    typename F,
	    typename B>
  void glbReduce(T* glb_res,
		 B& buf,
		 const int64_t nloc,
		 F&& f,
		 const int nslices=1,
		 const int nloc_slices=1,
		 const int loc_offset=0)
  {
    T loc_res[nslices];
    memset(loc_res,0,sizeof(T)*nslices);
    
    locReduce(loc_res+loc_offset,buf,nloc,nloc_slices,f);
    
    non_loc_reduce(glb_res,loc_res,nslices);
  }
}

// #undef EXTERN_REDUCE
// #undef INIT_REDUCE_TO

#endif
