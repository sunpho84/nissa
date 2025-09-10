#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/vectors.hpp>
#include <routines/ios.hpp>

#define EXTERN_BENCH
# include <base/bench.hpp>

#include <base/memory_manager.hpp>
#include <communicate/communicate.hpp>
#include <geometry/geometry_lx.hpp>

namespace nissa
{
  //copy memory
  void bench_memory_copy(double *out,double *in,int size)
  {
    size/=8;
    
    for(int i=0;i<size;i++) out[i]=in[i];
  }
  
  //benchmark memory
  void bench_memory_bandwidth(int mem_size)
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
    
    MASTER_PRINTF("time to copy %d Mbytes: %lg, %lg Mbs\n",mem_size/1024/1024,
		  bench_time,mem_size/1024.0/1024/bench_time);
  }
  
  //benchmark the net speed
  void benchNetSpeed()
  {
    MemoryManager *mem=memoryManager<defaultMemorySpace>();
    
    char* out=mem->provide<char>(sendBufSize);
    char* in=mem->provide<char>(sendBufSize);
    
    MASTER_PRINTF("Communication benchmark, packet size %lu\n",sendBufSize);
    fflush(stdout);
    
    //speeds
    const int ntests=10;
    
    for(int sRank=0;sRank<nranks;sRank++)
      for(int dRank=0;dRank<nranks;dRank++)
	if(sRank!=dRank)
	  if(sRank==rank or dRank==rank)
	    {
	      double speedAve=0,speed_var=0;
	      
	      for(int itest=0;itest<ntests;itest++)
		{
		  double time=-take_time();
		  const int tag=9;
		  MPI_Sendrecv(out,sendBufSize,MPI_CHAR,dRank,tag,in,sendBufSize,MPI_CHAR,sRank,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  time+=take_time();
		  
		  const double speed=sendBufSize/time/1e6;
		  
		  speedAve+=speed;
		  speed_var+=speed*speed;
		}
	      
	      //compute
	      speedAve/=ntests;
	      speed_var/=ntests;
	      speed_var-=speedAve*speedAve;
	      
	      const double speedStddev=
		sqrt(speed_var);
	      
	      printf("%d <---> %d : %lg, stddev %lg Mb/s\n",sRank,dRank,speedAve,speedStddev);
	      
	      MPI_Barrier(MPI_COMM_WORLD);
	    }
    
    mem->release(in);
    mem->release(out);
  }
}
