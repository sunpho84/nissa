#ifndef _UNROLL_HPP
#define _UNROLL_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#if defined(__NVCC__)||defined(__CLANG__)
# define UNROLL_PREFIX \
  _Pragma("unroll")
#else
# define UNROLL_PREFIX
#endif

#define UNROLL_FOR(I,MIN,MAX)			\
  UNROLL_PREFIX					\
  for(auto I=MIN;I<MAX;I++)

#define UNROLL_FOR_ALL_COLS(IC)			\
  UNROLL_FOR(IC,0,NCOL)

#define UNROLL_FOR_ALL_DIRS(MU)			\
  UNROLL_FOR(MU,0,NDIM)

#define UNROLL_FOR_ALL_DIRAC(ID)		\
  UNROLL_FOR(ID,0,NDIRAC)

namespace nissa
{
}

#endif
