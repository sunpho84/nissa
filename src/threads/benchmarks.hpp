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
    
    /// Creates from the name, file and line
    KernelInfo(const std::string& name) :
      name(name)
    {
    }
    
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
    KernelInfoLaunchParsStat(const std::string& name) :
      info(name),
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
    "autoTunedKernelInfo.snf";
  
  /// Stores the tuned kernel info
  inline void storeTunedKernelsInfo()
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
  inline void tryLoadTunedKernelsInfo()
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
    
    /// Read the number of tuned kernels
    int n;
    read_int(&n);
    master_printf("Found %d tuned kernels\n",n);
    
    for(int i=0;i<n;i++)
      {
	/// Read the kernel name
	char name[1024];
	read_str(name,1024);
	
	/// Read the number of tuned suzes
	int nSizes;
	read_int(&nSizes);
	
	kernelInfoLaunchParsStats.emplace_back(name);
	
	for(int i=0;i<nSizes;i++)
	  {
	    /// Read the loop length
	    int64_t loopLength;
	    read_int64_t(&loopLength);
	    
	    /// Read the optimal block size
	    int optimalBlockSize;
	    read_int(&optimalBlockSize);
	    
	    kernelInfoLaunchParsStats.back().launchParsStatList[loopLength].optimalBlockSize=optimalBlockSize;
	  }
      }
    
    verbosity_lv3_master_printf("-------- List of tuned kernels --------\n");
    int id=0;
    for(const auto& [kernelInfo,kernelLaunchParsStats] :  kernelInfoLaunchParsStats)
      {
	verbosity_lv3_master_printf("| %d - %s %zu\n",id,kernelInfo.name.c_str(),kernelLaunchParsStats.size());
	for([[maybe_unused]] const auto& u : kernelLaunchParsStats)
	  verbosity_lv3_master_printf("|| %ld %d\n",u.first,u.second.optimalBlockSize);
	id++;
      }
    verbosity_lv3_master_printf("-------------------\n");
    
    close_file(fin);
    input_global=back;
  }
  
  /// Reference to the flag specifying whether we are inside a parallel for
  extern int insideParallelFor;
  
  /// Fetches the optimal block size in the archive, or tunes the kernel
  template <typename F>
  int getOptimalBlockSize(const int& kernelId,
			  const int64_t& loopLength,
			  const int& minBlockSize,
			  const int& maxBlockSize,
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
	  printf("Benchmarking kernel %s for loop size %ld, doNotBackupDuringBenchmark: %d, insideParallelFor: %d\n",
		 info.name.c_str(),loopLength,doNotBackupDuringBenchmark,insideParallelFor);
	
	/// Takes note of the action to carry out before benchmarking
	std::vector<BenchmarkAction> benchmarkBeginActions(std::move(nissa::benchmarkBeginActions));
	
	/// Takes note of the action to carry out after benchmarking
	std::vector<BenchmarkAction> benchmarkEndActions(std::move(nissa::benchmarkEndActions));
    
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
	
	if(maxBlockSize<=minBlockSize)
	  crash("minbockSize %d cannot be equal or larger than maxBlockSize %d",minBlockSize,maxBlockSize);
	
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
	      (take_time()-initTime)/nBench;
	    
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
	
	storeTunedKernelsInfo();
      }
    else
      if(print)
	printf("Retrieved optimalBlockSize %d, nInvoke: %ld\n",
	       optimalBlockSize,launchParsStat.nInvoke);
    
    return optimalBlockSize;
  }
}

#endif
