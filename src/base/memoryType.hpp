#ifndef _MEMORYTYPE_HPP
#define _MEMORYTYPE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

namespace nissa
{
 enum class MemoryType{CPU=1,  ///< Memory allocated on CPU side
			GPU=2}; ///< Memory allocated on GPU side
  
  /// Default memory to be used
  constexpr MemoryType defaultMemoryType=
	      MemoryType::
#ifdef ENABLE_DEVICE_CODE
	      GPU
#else
	      CPU
#endif
	      ;
  
  /// GPU memory type if compiling for device, CPU otherwise
  constexpr MemoryType maybeGpuMemoryType=
	      MemoryType::
#ifdef ENABLE_DEVICE_CODE
	      GPU
#else
	      CPU
#endif
	      ;
}

#endif
