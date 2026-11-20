#ifndef _BENCHMARKS_HPP
#define _BENCHMARKS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cmath>
#include <cstdlib>
#include <functional>
#include <map>
#include <string>

#include "base/debug.hpp"
#include "io/input.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

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
    if(not fileExists(tunedKernelInfoFileName))
      {
	MASTER_PRINTF("%s file not found, ignoring load tuned parameters\n",tunedKernelInfoFileName);
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
    MASTER_PRINTF("Found %d tuned kernels\n",n);
    
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
    
    VERBOSITY_LV3_MASTER_PRINTF("-------- List of tuned kernels --------\n");
    int id=0;
    for(const auto& [kernelInfo,kernelLaunchParsStats] :  kernelInfoLaunchParsStats)
      {
	VERBOSITY_LV3_MASTER_PRINTF("| %d - %s %zu\n",id,kernelInfo.name.c_str(),kernelLaunchParsStats.size());
	for([[maybe_unused]] const auto& u : kernelLaunchParsStats)
	  VERBOSITY_LV3_MASTER_PRINTF("|| %ld %d\n",u.first,u.second.optimalBlockSize);
	id++;
      }
    VERBOSITY_LV3_MASTER_PRINTF("-------------------\n");
    
    close_file(fin);
    input_global=back;
  }
  
  /// Reference to the flag specifying whether we are inside a parallel for
  extern int insideParallelFor;
  
  /// Read the tuning trace switch once, sharing its cached value across kernel types
  inline bool getTraceLaunch()
  {
    static const bool traceLaunches=
      getExternalPar<bool>("NISSA_TUNING_TRACE","enable tracing of the kernel benchmark");
    
    return traceLaunches;
  }
  
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
    
    /// A reference block size
    const int defaultBlockSize=128;
    
    /// RunTime for default block size
    double defaultBlockRunTime=0.0;
    
    if(print)
      printf("ToBeTuned: %d\n",toBeTuned);
    
    if(toBeTuned)
      {
	/// Detailed launch diagnostics are opt-in because flushing perturbs tuning execution
	const bool traceLaunches=getTraceLaunch();
	
	/// Bracket each synchronous launch and measure its time without diagnostic output
	const auto measuredLaunch=
	  [&](const int blockSize,
	      const char* phase,
	      const int iteration)
	  {
	    if(traceLaunches)
	      {
		MASTER_PRINTF("NISSA tuning rank=%d kernel=\"%s\" loop=%ld block=%d grid=%ld phase=%s iteration=%d state=begin\n",
			      rank,info.name.c_str(),loopLength,blockSize,
			      loopLength>0?1+(loopLength-1)/blockSize:0,phase,iteration);
		fflush(stdout);
	      }
	    
	    const double start=take_time();
	    launch(blockSize);
	    const double elapsed=take_time()-start;
	    
	    if(traceLaunches)
	      {
		MASTER_PRINTF("NISSA tuning rank=%d kernel=\"%s\" loop=%ld block=%d grid=%ld phase=%s iteration=%d state=complete seconds=%.17g\n",
			      rank,info.name.c_str(),loopLength,blockSize,
			      loopLength>0?1+(loopLength-1)/blockSize:0,phase,iteration,elapsed);
		fflush(stdout);
	      }
	    
	    return elapsed;
	  };
	
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
	
	/// Repeat cheap candidates for precision and limit time spent on expensive ones
	const int minBench=3;
	const int maxBench=100;
	const double minMeasuredSeconds=0.2;
	const double maxMeasuredSeconds=2.0;
	const double maxRelativeStandardError=0.02;
	
	/// Minimal runtime for the kernel
	double minRunTime=0.0;
	
	if(maxBlockSize<=minBlockSize)
	  CRASH("minbockSize %d cannot be equal or larger than maxBlockSize %d",minBlockSize,maxBlockSize);
	
	if(print)
	  printf("starting test with block size of powers of two in the range [%d;%d)\n",minBlockSize,maxBlockSize);
	for(int testBlockSize=minBlockSize;testBlockSize<=maxBlockSize;testBlockSize*=2)
	  {
	    // warmup
	    measuredLaunch(testBlockSize,"warmup",0);
	    
	    /// Accumulate the launch times without storing them individually
	    RunningStatistics timing;
	    bool stable=false;
	    while(timing.nEstimates<maxBench and
		  (timing.nEstimates<minBench or (not stable and timing.sum<maxMeasuredSeconds)))
	      {
		const double elapsed=measuredLaunch(testBlockSize,"timed",timing.nEstimates);
		timing.push(elapsed);
		
		if(timing.nEstimates>=minBench and timing.sum>=minMeasuredSeconds and timing.average()>0.0)
		  stable=timing.standardError()<=maxRelativeStandardError*timing.average();
	      }
	    
	    const double runTime=timing.average();
	    
	    if(optimalBlockSize==0 or minRunTime>runTime)
	      {
		optimalBlockSize=testBlockSize;
		minRunTime=runTime;
	      }
	    
	    if(testBlockSize==defaultBlockSize)
	      defaultBlockRunTime=runTime;
	    
	    if(print)
	      printf("Benchmarked with blockSize %d, runtime %lg s from %zu launches (%.3g s total), minimal %lg s current optimal size %d\n",testBlockSize,runTime,timing.nEstimates,timing.sum,minRunTime,optimalBlockSize);
	  }
	
	if(not doNotBackupDuringBenchmark)
	  {
	    doNotBackupDuringBenchmark=true;
	    for(BenchmarkAction& e : benchmarkEndActions)
	      e();
	    doNotBackupDuringBenchmark=false;
	  }
	
	if(defaultBlockRunTime>0.0)
	  VERBOSITY_LV1_MASTER_PRINTF("Tuned kernel %s for loop length %ld, optimal block size: %d, minimal time: %lg s, time with default block size (%d) would have been %lg s, saved time: %lg %c\n",
				      info.name.c_str(),loopLength,optimalBlockSize,minRunTime,defaultBlockSize,defaultBlockRunTime,(1-minRunTime/defaultBlockRunTime)*100,'%');
	else
	  VERBOSITY_LV1_MASTER_PRINTF("Tuned kernel %s for loop length %ld, optimal block size: %d, minimal time: %lg s\n",
				      info.name.c_str(),loopLength,optimalBlockSize,minRunTime);
	
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
