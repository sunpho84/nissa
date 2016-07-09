#ifndef _INIT_HPP
#define _INIT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

namespace nissa
{
  const char compile_info[4][1024]={CONFIG_TIME,CONFIG_FLAGS,__TIME__,__DATE__};
  
  void init_nissa(int narg,char **arg,const char compile_info[5][1024]);
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024]);
  inline void init_nissa(int narg,char **arg)
  {init_nissa(narg,arg,compile_info);}
  inline void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg))
  {init_nissa_threaded(narg,arg,main_function,compile_info);}
}

#endif
