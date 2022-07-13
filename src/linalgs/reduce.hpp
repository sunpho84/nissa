#ifndef _REDUCE_HPP
#define _REDUCE_HPP

#include <cstdint>

#include "base/vectors.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#ifndef EXTERN_REDUCE
 #define EXTERN_REDUCE extern
 #define INIT_REDUCE_TO(...)
#else
 #define INIT_REDUCE_TO(...) __VA_ARGS__
#endif

namespace nissa
{
  /// Size of the reducing buffer
  EXTERN_REDUCE int64_t reducing_buffer_size INIT_REDUCE_TO(=0);
  
  /// Reference to the reducing buffer
  EXTERN_REDUCE CUDA_MANAGED void *reducing_buffer INIT_REDUCE_TO(=nullptr);
  
  /// Get the reduction buffer
  template <typename T>
  T* get_reducing_buffer(const int64_t& n)
  {
    const int64_t required_size=sizeof(T)*n;
    if(reducing_buffer_size<required_size)
      {
	if(reducing_buffer!=nullptr)
	  nissa_free(reducing_buffer);
	
	//reducing_buffer_size=0; check
      }
    
    if(reducing_buffer==nullptr)
      {
	if(reducing_buffer_size==0)
	  master_printf("Allocating the reduction buffer to %ld bytes\n",required_size);
	else
	  master_printf("Reallocating the reduction buffer from %ld bytes to %ld bytes\n",reducing_buffer_size,required_size);
	
	reducing_buffer=nissa_malloc("reducing_buffer",required_size,char);
	reducing_buffer_size=required_size;
      }
    
    return (T*)reducing_buffer;
  }
  
  void deallocate_reduction_buffer();
  
  void loc_reduce(int64_t* loc_res,int64_t* buf,int64_t n,int nslices=1);
  void loc_reduce(double* loc_res,double* buf,int64_t n,int nslices=1);
  void loc_reduce(complex* loc_res,complex* buf,int64_t n,int nslices=1);
#ifdef REPRODUCIBLE_RUN
  void loc_reduce(float_128* loc_res,float_128* buf,int64_t n,int nslices=1);
  void loc_reduce(complex_128* loc_res,complex_128* buf,int64_t n,int nslices=1);
#endif
  
  //reduce a vector
  template <typename T>
  inline void non_loc_reduce(T* out_glb,T* in_loc=nullptr,const int n=1,MPI_Op mpi_op=MPI_Op_sum_for_type<T>())
  {
    if(in_loc==nullptr)
      in_loc=(T*)MPI_IN_PLACE;
    
    MPI_Allreduce(in_loc,out_glb,n,MPI_Datatype_of<T>(),mpi_op,MPI_COMM_WORLD);
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T>
  CUDA_HOST_AND_DEVICE
  void reduceSummer(T& out,const T& in)
  {
    out+=in;
  }
  
  template <>
  CUDA_HOST_AND_DEVICE
  inline void reduceSummer(complex& out,const complex& in)
  {
    complex_summassign(out,in);
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T>
  CUDA_HOST_AND_DEVICE
  void reduceAssigner(T& out,const T& in)
  {
    out=in;
  }
  
  template <>
  CUDA_HOST_AND_DEVICE
  inline void reduceAssigner(complex& out,const complex& in)
  {
    complex_copy(out,in);
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T,
	    typename F>
  void locReduce(T *loc_res,T *buf,int64_t n,const int nslices,F&& f)
  {
    master_printf("%s function\n",__PRETTY_FUNCTION__);
    
    if(n%nslices)
      crash("number of elements %ld not divisible by number of slices %d",n,nslices);
    
    const int64_t nori_per_slice=n/nslices;
    int64_t nper_slice=n/nslices;
    verbosity_lv2_master_printf("n: %lld, nslices: %d, nori_per_slice: %lld nper_slice: %ld\n",nslices,nori_per_slice,nper_slice);
    
    const double init_time=take_time();
    while(nper_slice>1)
      {
	const int64_t stride=(nper_slice+1)/2;
	const int64_t nreductions_per_slice=nper_slice/2;
	const int64_t nreductions=nreductions_per_slice*nslices;
	verbosity_lv3_master_printf("nper_slice: %lld, stride: %lld, nreductions_per_slice: %lld, nreductions: %lld \n",nper_slice,stride,nreductions_per_slice,nreductions);
	
	
	NISSA_PARALLEL_LOOP(ireduction,0,nreductions)
	  {
	    const int64_t islice=ireduction%nslices;
	    const int64_t ireduction_in_slice=ireduction/nslices;
	    const int64_t first=ireduction_in_slice+nori_per_slice*islice;
	    const int64_t second=first+stride;
	    
	    f(buf[first],buf[second]);
	  }
	  NISSA_PARALLEL_LOOP_END;
	  THREAD_BARRIER();
	  
	  nper_slice=stride;
      }
    
    master_printf("reduction ended, took %lg s\n",take_time()-init_time);
    
    for(int islice=0;islice<nslices;islice++)
      reduceAssigner(loc_res[islice],buf[islice*nori_per_slice]);
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Reduce a vector over all nodes
  template <typename T>
  void glb_reduce(T* glb_res,T* buf,int64_t nloc,const int nslices=1,const int nloc_slices=1,const int loc_offset=0)
  {
    T loc_res[nslices];
    memset(loc_res,0,sizeof(T)*nslices);
    
    locReduce(loc_res+loc_offset,buf,nloc,nloc_slices,[] CUDA_DEVICE (auto& res,const auto& acc)  __attribute__((always_inline))
    {
      reduceSummer(res,acc);
    });
    
    non_loc_reduce(glb_res,loc_res,nslices);
  }
}

#undef EXTERN_REDUCE
#undef INIT_REDUCE_TO

#endif
