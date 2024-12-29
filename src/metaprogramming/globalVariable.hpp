#ifndef _GLOBALVARIABLE_HPP
#define _GLOBALVARIABLE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <metaprogramming/inline.hpp>

namespace nissa
{
#if defined(USE_CUDA) && defined(__CUDA_ARCH__)
# define PROVIDE_CONST_ACCESSOR_FOR_GLOBAL_VAR(TYPE,NAME)	\
  inline __device__ const TYPE& NAME=hidden::_## NAME ## Gpu
#else
# define PROVIDE_CONST_ACCESSOR_FOR_GLOBAL_VAR(TYPE,NAME)	\
  inline const TYPE& NAME=hidden::_## NAME ## Cpu
#endif

#ifdef USE_CUDA
# define MAYBE_PROVIDE_DEVICE_GLOBAL_VAR(TYPE,NAME)		\
  __device__ __constant__ inline TYPE _ ## NAME ## Gpu
# define MAYBE_COPY_GLOBAL_VAR_TO_DEVICE(TYPE,NAME)		\
  cudaMemcpyToSymbol(hidden::_ ## NAME ## Gpu,			\
		     &hidden::_ ## NAME ## Cpu,			\
		     sizeof(TYPE))

#else
# define MAYBE_PROVIDE_DEVICE_GLOBAL_VAR(TYPE,NAME)
# define MAYBE_COPY_GLOBAL_VAR_TO_DEVICE(TYPE,NAME)
#endif

#define PROVIDE_GLOBAL_VAR_STORAGE(TYPE,NAME)		\
  namespace hidden					\
  {							\
    inline TYPE _ ## NAME ## Cpu;			\
							\
    MAYBE_PROVIDE_DEVICE_GLOBAL_VAR(TYPE,NAME);		\
  }

#define PROVIDE_GLOBAL_VAR_SETTER(TYPE,NAME)		\
  INLINE_FUNCTION void set_ ## NAME(const TYPE& val)	\
  {						\
    using namespace hidden;			\
						\
    _ ## NAME ## Cpu=val;			\
    						\
    MAYBE_COPY_GLOBAL_VAR_TO_DEVICE(TYPE,NAME);	\
  }

#define PROVIDE_GLOBAL_VAR(TYPE,NAME)			\
  PROVIDE_GLOBAL_VAR_STORAGE(TYPE,NAME);		\
  PROVIDE_CONST_ACCESSOR_FOR_GLOBAL_VAR(TYPE,NAME);	\
  PROVIDE_GLOBAL_VAR_SETTER(TYPE,NAME)
}

#endif
