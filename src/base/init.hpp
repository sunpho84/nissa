#ifndef _INIT_HPP
#define _INIT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

namespace nissa
{
  const char compile_info[4][1024]={CONFIG_TIME,CONFIG_FLAGS,__TIME__,__DATE__};
  
  int bulk_recip_lat_volume(int *P,int *L);
  int bulk_volume(int *L);
  int compute_border_variance(int *L,int *X,int consider_reciprocal);
  void find_minimal_surface_grid(int *mP,int *L,int NP);
  void init_grid(int T,int L);
  void init_nissa(int narg,char **arg,const char compile_info[5][1024]);
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024]);
  inline void init_nissa(int narg,char **arg)
  {init_nissa(narg,arg,compile_info);}
  inline void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg))
  {init_nissa_threaded(narg,arg,main_function,compile_info);}
}

#endif
