#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <chrono>

#define EXTERN_BENCH
# include <base/bench.hpp>

#include <base/field.hpp>
#include <base/vectors.hpp>
#include <base/memory_manager.hpp>
#include <communicate/communicate.hpp>
#include <geometry/geometry_lx.hpp>
#include <new_types/su3_op.hpp>
#include <routines/ios.hpp>

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
  
  template <typename F>
  void  bench(F&& f,
	      const char* what,
	      const size_t& nTests)
  {
    /// Benchmark statistics
    struct Stat
    {
      double sum{};
      
      double sum2{};
      
      int nTests{};
      
      Stat& operator+=(const double& t)
      {
	sum+=t;
	sum2+=t*t;
	nTests++;
	
	return *this;
      }
      
      std::pair<double,double> get() const
      {
	const double ave=sum/nTests;
	const double var=sum2/nTests-ave*ave;
	
	return {ave,sqrt(var/(nTests-1))};
      }
    };
    
    Stat glb;
    
    Stat loc;
    
    //warmup
    f();
    
    for(size_t i=0;i<nTests;i++)
      {
	const double initMoment=take_time();
	std::chrono::time_point locInitMoment=std::chrono::steady_clock::now();
        f();
	auto locT=locInitMoment-std::chrono::steady_clock::now();
	const double t=take_time()-initMoment;
	
	loc+=std::chrono::duration_cast<std::chrono::nanoseconds>(locT).count()/1e9;
	glb+=t;
      }
    
    auto g=glb.get();
    
    auto l=loc.get();
    
    MASTER_PRINTF("Time to %s: %lg s +/- %lg s\n",
		  what,g.first,g.second);
    for(size_t iRank=0;iRank<nranks;iRank++)
      {
	if(rank==iRank)
	  {
	    printf(" Local estimate on rank %zu: %lg s +/- %lg s\n",
		   iRank,l.first,l.second);
	    fflush(stdout);
	  }
	MPI_Barrier(MPI_COMM_WORLD);
      }
  }
  
  /// Benchmark halo exchange
  void benchHaloExchange()
  {
    LxField<spincolor> f("f",WITH_HALO);
    f.reset();
    
    const size_t nTests=8;
    bench([&f]()
    {
      f.updateHalo();
    },"update halo on a spincolor lx field",nTests);
  }
  
  /// Benchmark the covariant shift
  void benchCovariantShift()
  {
    LxField<spincolor> out("out");
    LxField<quad_su3> conf("conf");
    LxField<spincolor> in("in",WITH_HALO);
    
    conf.reset();
    in.reset();
    
    in.updateHalo();
    
    const size_t nTests=8;
    bench([&in,
	   &conf,
	   &out]()
    {
      PAR(0,
	  locVol,
	  CAPTURE(TO_READ(in),
		  TO_READ(conf),
		  TO_WRITE(out)),
	  ivol,
	  {
	    unsafe_su3_prod_spincolor(out[ivol],conf[ivol][0],in[loclxNeighup[ivol][0]]);
	  });
    },
	  "perform a covariant shift of a spincolor field",
	  nTests);
  }
}
