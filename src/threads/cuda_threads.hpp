#ifndef _CUDA_THREADS_HPP
#define _CUDA_THREADS_HPP

#include "base/init.hpp"

#define NUM_THREADS 128

#define NACTIVE_THREADS 1
#define MANDATORY_PARALLEL
#define MANDATORY_NOT_PARALLEL
#define IS_PARALLEL false
#define GET_THREAD_ID()
#define THREAD_ID (0)
#define THREAD_BARRIER_FORCE()
#define THREAD_BARRIER()
#define IS_MASTER_THREAD (1)
#define NISSA_PARALLEL_LOOP(INDEX,EXT_START,EXT_END) cuda_parallel_for(EXT_START,EXT_END,[=] __host__ __device__ (const uint64_t& INDEX){
#define NISSA_PARALLEL_LOOP_END })
#define THREAD_ATOMIC_EXEC(inst) inst
#define THREAD_BROADCAST(out,in) (out)=(in)
#define THREAD_BROADCAST_PTR(out,in) THREAD_BROADCAST(out,in)
#define THREADABLE_FUNCTION_0ARG(FUNC_NAME) void FUNC_NAME(){
#define THREADABLE_FUNCTION_1ARG(FUNC_NAME,AT1,A1) void FUNC_NAME(AT1 A1){
#define THREADABLE_FUNCTION_2ARG(FUNC_NAME,AT1,A1,AT2,A2) void FUNC_NAME(AT1 A1,AT2 A2){
#define THREADABLE_FUNCTION_3ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3){
#define THREADABLE_FUNCTION_4ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4){
#define THREADABLE_FUNCTION_5ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5){
#define THREADABLE_FUNCTION_6ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6){
#define THREADABLE_FUNCTION_7ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7){
#define THREADABLE_FUNCTION_8ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8){
#define THREADABLE_FUNCTION_9ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9){
#define THREADABLE_FUNCTION_10ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10){
#define THREADABLE_FUNCTION_11ARG(FUNC_NAME,AT1,A1,AT2,A2,AT3,A3,AT4,A4,AT5,A5,AT6,A6,AT7,A7,AT8,A8,AT9,A9,AT10,A10,AT11,A11) void FUNC_NAME(AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,AT6 A6,AT7 A7,AT8 A8,AT9 A9,AT10 A10,AT11 A11){

namespace nissa
{
  template <typename IMin,
	    typename IMax,
	    typename F>
  __global__
  void cuda_generic_kernel(const IMin &min,
			   const IMax &max,
			   F f)
  {
    const auto i=min+blockIdx.x*blockDim.x+threadIdx.x;
    if(i<max)
      f(i);
  }
  
  template <typename IMin,
	    typename IMax,
	    typename F>
  void cuda_parallel_for(const IMin &min,
			 const IMax &max,
			 F f)
  {
    const auto length=(max-min);
    const dim3 block_dimension(NUM_THREADS);
    const dim3 grid_dimension((length+block_dimension.x-1)/block_dimension.x);
    cuda_generic_kernel<<<grid_dimension,block_dimension>>>(min,max,f);
  }
  
  inline void cache_flush()
  {
  }
  
  inline void thread_barrier_internal()
  {
    cudaDeviceSynchronize();
  }
  
  //start nissa
  inline void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024])
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg,compile_info);
    
    main_function(narg,arg);
  }
}

#endif
