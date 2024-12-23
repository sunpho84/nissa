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
  inline void storeKernelInfo()
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
  inline void tryLoadKernelInfo()
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
}

#endif
