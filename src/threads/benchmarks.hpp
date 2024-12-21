#ifndef _BENCHMARKS_HPP
#define _BENCHMARKS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <functional>
#include <map>
#include <string>

namespace nissa
{
  using BenchmarkAction=
    std::function<void(void)>;
  
  using BenchmarkActions=
    std::vector<std::function<void(void)>>;
  
  inline BenchmarkActions benchmarkBeginActions;
  
  inline BenchmarkActions benchmarkEndActions;
  
  inline bool benchmarkInProgress;
  
  struct KernelInfo
  {
    std::string name;
    
    std::string file;
    
    int line;
    
    KernelInfo(const std::string& name,
	       const std::string& file,
	       const int& line) :
      name(name),
      file(file),
      line(line)
    {
    }
    
    KernelInfo()=default;
    
    KernelInfo(KernelInfo&&)=default;
  };
  
  struct KernelSizeLaunchParsStat
  {
    double totalTime{0.0};
    
    size_t nInvoke{0};
    
    int optimalBlockSize{0};
  };
  
  struct KernelInfoLaunchParsStat
  {
    KernelInfo info;
    
    std::map<int64_t,KernelSizeLaunchParsStat> launchParsStatList;
    
    KernelInfoLaunchParsStat(const std::string& name,
			     const std::string& file,
			     const int& line) :
      info(name,file,line),
      launchParsStatList()
    {
    }
  };
  
  inline std::vector<KernelInfoLaunchParsStat> kernelInfoLaunchParsStats;
}

#endif
