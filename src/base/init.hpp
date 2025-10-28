#ifndef _INIT_HPP
#define _INIT_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <geometry/geometry_lx.hpp>

namespace nissa
{
  /// Info on compilation and configuration
  const char compileConfigInfo[4][1024]=
    {CONFIG_TIME,CONFIG_FLAGS,__TIME__,__DATE__};
  
  int bulk_recip_lat_volume(int *P,int *L);
  int bulk_volume(int *L);
  int compute_border_variance(int *L,int *X,int consider_reciprocal);
  
  Coords findMinimalSurfaceGrid(const Coords& ext_L,
				const int& NR);
  
  void initGrid(const int& T,
		 const int& L);
  
  void initNissa(int narg,
		 char **arg,
		 const char compileConfigInfo[5][1024]);
  
  inline bool nissaInited;
  
  inline void initNissa(int narg,
			char **arg)
  {
    initNissa(narg,arg,compileConfigInfo);
  }
}

#endif
