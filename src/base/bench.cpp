#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "new_types/float_128.hpp"
#include "routines/ios.hpp"

#include "threads/threads.hpp"

#define EXTERN_BENCH
 #include "bench.hpp"

#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"

namespace nissa
{
  //copy memory
  void bench_memory_copy(double *out,double *in,int size)
  {
    size/=8;
    
    NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,size,THREAD_ID,NACTIVE_THREADS);
    
    for(int i=start;i<end;i++) out[i]=in[i];
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
    
    master_printf("time to copy %d Mbytes: %lg, %lg Mbs\n",mem_size/1024/1024,
		  bench_time,mem_size/1024.0/1024/bench_time);
  }
  
  //benchmark the net speed
  void bench_net_speed()
  {
    if(nranks>1)
      for(int ipow=14;ipow<=log2(send_buf_size);ipow+=2)
	{
	  //allocate a buffer
	  int size=1<<ipow;
	  char *out,*in;
	  
	  out=nissa_malloc("out",size,char);
	  in=nissa_malloc("in",size,char);
	  
	  //total time
	  double tot_time=0;
	  //speeds
	  double speed_ave=0,speed_var=0;
	  int ntests=10,n=0;
	  for(int itest=0;itest<ntests;itest++)
	    for(int drank=1;drank<nranks;drank++)
	      {
		double time=-take_time();
		int tag=9;
		MPI_Sendrecv(out,size,MPI_CHAR,(rank+drank)%nranks,tag,in,size,MPI_CHAR,(rank-drank+nranks)%nranks,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		time+=take_time();
		double speed=size/time/1e6;
		
		//increase
		n++;
		speed_ave+=speed;
		speed_var+=speed*speed;
		
		tot_time+=time;
	      }
	  
	  //reduce
	  MPI_Allreduce(MPI_IN_PLACE,&n,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	  MPI_Allreduce(MPI_IN_PLACE,&speed_ave,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  MPI_Allreduce(MPI_IN_PLACE,&speed_var,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  
	  //compute
	  speed_ave/=n;
	  speed_var/=n;
	  speed_var-=speed_ave*speed_ave;
	  double speed_stddev=sqrt(speed_var);
	  
	  nissa_free(in);
	  nissa_free(out);
	  
	  master_printf("Communication benchmark, packet size %d (%lg, stddev %lg) Mb/s (%lg s total)\n",size,speed_ave,speed_stddev,tot_time);
	}
  }
}
