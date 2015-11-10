#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //separate the even and odd part of a vector
  THREADABLE_FUNCTION_3ARG(split_lx_vector_into_eo_parts_internal, char**,out_eo, char*,in_lx, size_t,bps)
  {
    GET_THREAD_ID();
    
#ifdef BENCH
    if(IS_MASTER_THREAD)
      {
	remap_time-=take_time();
	nremap++;
      }
#endif
    
    //split
    NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
      memcpy(out_eo[loclx_parity[loclx]]+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
    
#ifdef BENCH
    if(IS_MASTER_THREAD) remap_time+=take_time();
#endif
    
    set_borders_invalid(out_eo[0]);
    set_borders_invalid(out_eo[1]);
  }
  THREADABLE_FUNCTION_END
  
  //paste the even and odd parts of a vector into a full lx vector
  THREADABLE_FUNCTION_3ARG(paste_eo_parts_into_lx_vector_internal, char*,out_lx, char**,in_eo, size_t,bps)
  {
    GET_THREAD_ID();
    
#ifdef BENCH
    if(IS_MASTER_THREAD)
      {
	remap_time-=take_time();
	nremap++;
      }
#endif
    
    //paste
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(eo,0,loc_volh)
	memcpy(out_lx+bps*loclx_of_loceo[par][eo],in_eo[par]+bps*eo,bps);
    
#ifdef BENCH
    if(IS_MASTER_THREAD) remap_time+=take_time();
#endif
    
    set_borders_invalid(out_lx);
  }
  THREADABLE_FUNCTION_END
}
