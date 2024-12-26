#ifndef _BENCHMARKS_HPP
#define _BENCHMARKS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cmath>
#include <functional>
#include <map>
#include <string>

#include "base/debug.hpp"
#include "io/input.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  /// Type representing an action to be carried out at benchmark begin/end
  using BenchmarkAction=
    std::function<void(void)>;
  
  /// Type representing the list of actions to be carried out at benchmark begin/end
  using BenchmarkActions=
    std::vector<std::function<void(void)>>;
  
  /// List of actions to be carried out at benchmark begin
  inline BenchmarkActions benchmarkBeginActions;
  
  /// List of actions to be carried out at benchmark end
  inline BenchmarkActions benchmarkEndActions;
  
  /// States whether needs to backup before doing the benchmark
  ///
  // Needed to bypass nested backup
  inline bool doNotBackupDuringBenchmark;
  
  /// Stores the info of a kernel
  struct KernelInfo
  {
    /// Name of the kernel
    std::string name;
    
    /// File where the kernel is defined
    std::string file;
    
    /// Line where the kernel is defined
    int line;
    
    /// Creates from the name, file and line
    KernelInfo(const std::string& name,
	       const std::string& file,
	       const int& line) :
      name(name),
      file(file),
      line(line)
    {
    }
    
    /// Default constructor
    KernelInfo()=default;
    
    /// Default move constructor
    KernelInfo(KernelInfo&&)=default;
  };
  
  /// Stores the statistics of execution of a kernel, and the optimal block size to execute
  struct KernelSizeLaunchParsStat
  {
    /// Total execution time
    double totalTime{0.0};
    
    /// Number of invokation
    size_t nInvoke{0};
    
    /// Optimal block size for execution
    int optimalBlockSize{0};
  };
  
  /// Info, statistics and parameters to launch a kernel
  struct KernelInfoLaunchParsStat
  {
    /// Kernel info
    KernelInfo info;
    
    /// Statistics of execution of each individiual launch size
    std::map<int64_t,KernelSizeLaunchParsStat> launchParsStatList;
    
    /// Constructor from name, file and line
    KernelInfoLaunchParsStat(const std::string& name,
			     const std::string& file,
			     const int& line) :
      info(name,file,line),
      launchParsStatList()
    {
    }
  };
  
  /// Stores the info, statistics and parameters to launch a kernel
  inline std::vector<KernelInfoLaunchParsStat> kernelInfoLaunchParsStats;
  
  /// Overhead of the benchmark
  inline std::pair<double,double> benchOverhead;
  
  /// Compute and stores the benchmark overhead
  inline std::pair<double,double> estimateBenchOverhead()
  {
    double x=0.0;
    
    double x2=0.0;
    
    const int nt=1000;
    for(int t=0;t<nt;t++)
      {
	const double i=take_time();
	
	const double e=take_time();
	
	const double d=e-i;
	
	x+=d;
	
	x2+=d*d;
      }
    
    x/=nt;
    
    x2/=nt;
    
    x2-=x*x;
    
    return std::make_pair(x,sqrt(x2/(nt-1)));
  }
  
  /// Name of the file used to store the tuned kernel info
  constexpr char tunedKernelInfoFileName[]=
    "tunedKernelInfo.txt";
  
  /// Stores the tuned kernel info
  inline void storeKernelTunedInfo()
  {
    FILE* fout=open_file(tunedKernelInfoFileName,"w");
    
    master_fprintf(fout,"%zu\n",kernelInfoLaunchParsStats.size());
    for(const auto& [kernelInfo,kernelLaunchParsStats] :  kernelInfoLaunchParsStats)
      {
	master_fprintf(fout,"%s %zu\n",kernelInfo.name.c_str(),kernelLaunchParsStats.size());
	for(const auto& u : kernelLaunchParsStats)
	  master_fprintf(fout,"%ld %d\n",u.first,u.second.optimalBlockSize);
      }
    
    close_file(fout);
  }
  
  /// Tries to load the tuned kernel info
  inline void tryLoadKernelTunedInfo()
  {
    if(not file_exists(tunedKernelInfoFileName))
      {
	master_printf("%s file not found, ignoring load tuned parameters\n",tunedKernelInfoFileName);
	return ;
      }
    
    FILE* fin=
      open_file(tunedKernelInfoFileName,"r");
    
    /// hack
    FILE* back=input_global;
    input_global=fin;
    
    int n;
    read_int(&n);
    master_printf("Found %d tuned kernels\n",n);
    
    for(int i=0;i<n;i++)
      {
	char name[1024];
	read_str(name,1024);
	
	int nSizes;
	read_int(&nSizes);
	
	kernelInfoLaunchParsStats.emplace_back(name,"",0);
	
	for(int i=0;i<nSizes;i++)
	  {
	    int64_t size;
	    read_int64_t(&size);
	    
	    int optimalBlockSize;
	    read_int(&optimalBlockSize);
	    
	    kernelInfoLaunchParsStats.back().launchParsStatList[size].optimalBlockSize=optimalBlockSize;
	  }
      }
    
    close_file(fin);
    input_global=back;
  }
  
  /// Fetches the optimal block size in the archive, or tunes the kernel
  template <typename F>
  int getOptimalBlockSize(const int& kernelId,
			  const int64_t& loopLength,
			  const int minBlockSize,
			  const int maxBlockSize,
			  const F& launch)
  {
    /// Decide whether to print
    const bool print=
      (VERBOSITY_LV3
       and rank==0);
    
    /// Fetch the info and the parameters of the kernel launches
    auto& [info,launchParsStatList]=
      kernelInfoLaunchParsStats[kernelId];
    
    /// Retrieves the parameters to launch wih thte given loop size
    const auto re=
      launchParsStatList.try_emplace(loopLength);
    
    /// Determine if tuned needs to be done or not
    bool toBeTuned=
      re.second;
    
    /// Reference to the parameters for the launch
    KernelSizeLaunchParsStat& launchParsStat=
      re.first->second;
    
    /// Reference to the optimal block size
    int& optimalBlockSize=
      launchParsStat.optimalBlockSize;
    
    if(0)
      {
	optimalBlockSize=128;
	toBeTuned=false;
      }
    
    if(print)
      printf("ToBeTuned: %d\n",toBeTuned);
    
    if(toBeTuned)
      {
	if(print)
	  printf("Benchmarking kernel %s for loop size %ld\n",
		 info.name.c_str(),loopLength);
	
	if(not doNotBackupDuringBenchmark)
	  {
	    doNotBackupDuringBenchmark=true;
	    for(BenchmarkAction& b : benchmarkBeginActions)
	      b();
	    doNotBackupDuringBenchmark=false;
	  }
	
	optimalBlockSize=0;
	
	const int nBench=100;
	double minTime=0.0;
	
	if(print)
	  printf("starting test with block size of powers of two in the range [%d;%d)\n",minBlockSize,maxBlockSize);
	for(int testBlockSize=minBlockSize;testBlockSize<=maxBlockSize;testBlockSize*=2)
	  {
	    // warmup
	    launch(testBlockSize);
	    
	    /// Initial time
	    const double initTime=
	      take_time();
	    
	    for(int i=0;i<nBench;i++)
	      launch(testBlockSize);
	    
	    /// Execution time
	    const double runTime=
	      take_time()-initTime;
	    
	    if(optimalBlockSize==0 or minTime>runTime)
	      {
		optimalBlockSize=testBlockSize;
		minTime=runTime;
	      }
	    
	    if(print)
	      printf("Benchmarked with blockSize %d, runtime %lg s minimal %lg s current optimal size %d\n",testBlockSize,runTime,minTime,optimalBlockSize);
	  }
	
	if(not doNotBackupDuringBenchmark)
	  {
	    doNotBackupDuringBenchmark=true;
	    for(BenchmarkAction& e : benchmarkEndActions)
	      e();
	    doNotBackupDuringBenchmark=false;
	  }
	
	storeKernelTunedInfo();
      }
    else
      if(print)
	printf("Retrieved optimalBlockSize %d, nInvoke: %ld\n",
	       optimalBlockSize,launchParsStat.nInvoke);
    
    return optimalBlockSize;
  }
}

#endif
