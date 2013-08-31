#include "../routines/ios.h"
#ifdef USE_THREADS
 #include "../routines/thread.h"
#endif
#ifdef BGQ
 #include "../bgq/intrinsic.h"
#endif

#include "global_variables.cpp"
#include "thread_macros.h"

//copy memory
void bench_memory_copy(double *out,double *in,int size)
{
  GET_THREAD_ID();
  size/=8;

  NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,size,thread_id,NACTIVE_THREADS);

#if  BGQ
  double *temp_out=out-4;
  double *temp_in=in-4;
  for(int i=start;i<end;i+=4)
    {
      reg_bi_complex reg;
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(reg,temp_in);
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(temp_out,reg);
    }
#else
  for(int i=start;i<end;i++) out[i]=in[i];
#endif
}

//benchmark memory
THREADABLE_FUNCTION_1ARG(bench_memory_bandwidth, int,mem_size)
{
  //allocate double
  double *a=nissa_malloc("a",mem_size/sizeof(double),double);
  double *b=nissa_malloc("b",mem_size/sizeof(double),double);
  
  //first call to warm up
  bench_memory_copy(a,b,mem_size);

  //exec 10 times
  int ntests=10;
  double bench_time=-take_time();
  for(int i=0;i<ntests;i++) bench_memory_copy(a,b,mem_size);
  bench_time+=take_time();
  bench_time/=ntests;

  nissa_free(a);
  nissa_free(b);

  master_printf("time to copy %d Mbytes: %lg, %lg Mbs\n",mem_size/1024/1024,bench_time,mem_size/1024/1024/bench_time);
}}
