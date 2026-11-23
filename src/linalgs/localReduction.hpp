#ifndef _LOCAL_REDUCTION_HPP
#define _LOCAL_REDUCTION_HPP

#include <charconv>
#include <cstdlib>
#include <string_view>

#include "localReductionV1.hpp"
#include "localReductionV2.hpp"

namespace nissa
{
  /// Local algorithm only; the global MPI collective is shared by both versions
  enum class LocalReductionVersion{V1=1,V2=2};
  
  /// Runtime selector, allowing comparison with the same executable and input
  inline constexpr char localReductionVersionEnvironment[]="NISSA_LOCAL_REDUCTION_VERSION";
  
  /// Read once per process, keeping the original implementation as the default
  inline LocalReductionVersion localReductionVersion()
  {
    static const auto version=[]()
    {
      int selected=1;
      getExternalPar(selected,localReductionVersionEnvironment,"Version for local reduction, 1 (default) or 2");
      if(selected!=1 and selected!=2)
	CRASH("%s must be 1 or 2",localReductionVersionEnvironment);
      
      return static_cast<LocalReductionVersion>(selected);
    }();
    
    return version;
  }
  
  /// Reduce equal nonempty slices using the selected associative, commutative operation
  /// V1 consumes its input; V2 preserves it except for its explicitly reported type fallback
  template <typename T,
            typename B,
            typename Op>
  void locReduce(T* result,
                 B& buffer,
                 const int64_t n,
                 const int nslices,
                 Op&& op,
                 const LocalReductionVersion version=localReductionVersion())
  {
    if(nslices<=0 or n<=0 or n%nslices!=0)
      CRASH("Local reduction requires positive n and slices with n divisible by slices (n=%ld, slices=%d)",n,nslices);
    
    if(version==LocalReductionVersion::V1)
      locReduceV1(result,buffer,n,nslices,op);
    else if(version==LocalReductionVersion::V2)
      locReduceV2(result,buffer,n,nslices,op);
    else
      CRASH("Unsupported local reduction version");
  }
}

#endif
