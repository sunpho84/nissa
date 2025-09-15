#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
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
    const size_t size=sendBufSize/10;
    
    MemoryManager *mem=memoryManager<defaultMemorySpace>();
    
    char* out=mem->provide<char>(size);
    char* in=mem->provide<char>(size);
    
    MASTER_PRINTF("Communication benchmark, packet size: %lu bytes\n",size);
    fflush(stdout);
    
    //speeds
    const int ntests=4;
    
    for(int sRank=0;sRank<nranks;sRank+=nloc_ranks)
      for(int dRank=0;dRank<nranks;dRank+=nloc_ranks)
	if(sRank!=dRank)
	  {
	    if(sRank==rank or dRank==rank)
	      {
		double speedAve=0,speedVar=0;
		
		for(int itest=0;itest<ntests;itest++)
		  {
		    double time=-take_time();
		    const int tag=9;
		    if(rank==sRank)
		      {
			// printf("on rank %d going to send %lu to %d\n",rank,size,dRank);
			// fflush(stdout);
			MPI_Send(out,size,MPI_CHAR,dRank,tag,MPI_COMM_WORLD);
		      }
		    else
		      {
			// printf("on rank %d going to receive %lu from %d\n",rank,size,sRank);
			// fflush(stdout);
			MPI_Recv(in,size,MPI_CHAR,sRank,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		      }
		    
		    time+=take_time();
		    
		    const double speed=size/time/1e9;
		    
		    speedAve+=speed;
		    speedVar+=speed*speed;
		  }
		
		//compute
		speedAve/=ntests;
		speedVar/=ntests;
		speedVar-=speedAve*speedAve;
		
		const double speedStddev=
		  sqrt(speedVar);
		
		if(sRank==rank)
		  printf("%d ---> %d : %lg, stddev %lg GB/s\n",sRank,dRank,speedAve,speedStddev);
		fflush(stdout);
	      }
	    
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
    
    mem->release(in);
    mem->release(out);
  }
  
  /// Benchmark halo exchange
  void benchHaloExchange()
  {
    LxField<spincolor> f("f",WITH_HALO);
    f.reset();
    
    double timeAve{};
    double timeVar{};
    
    const size_t nTests=8;
    for(size_t i=0;i<nTests;i++)
      {
	const double initMoment=take_time();
	f.updateHalo(true);
	const double t=take_time()-initMoment;
	timeAve+=t;
	timeVar+=t*t;
      }
    timeAve/=nTests;
    timeVar/=nTests;
    timeVar-=timeAve*timeAve;
    
    MASTER_PRINTF("Time to update halo on a spincolor lx field: %lg s +/- %lg s\n",timeAve,timeVar);
  }
}
