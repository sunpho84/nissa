#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/bench.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry_eo.hpp"
#include "geometry_lx.hpp"
#include "geometry_Leb.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //separate the even and odd part of a vector
  THREADABLE_FUNCTION_3ARG(split_lx_vector_into_eo_parts_internal, char**,out_eo, char*,in_lx, size_t,bps)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    //split
    NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
      memcpy(out_eo[loclx_parity[loclx]]+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out_eo[0]);
    set_borders_invalid(out_eo[1]);
  }
  THREADABLE_FUNCTION_END
  
  //separate the even and odd part of a vector
  THREADABLE_FUNCTION_4ARG(get_evn_or_odd_part_of_lx_vector_internal, char*,out_ev_or_od, char*,in_lx, size_t,bps, int,par)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    //get
    NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
      if(loclx_parity[loclx]==par)
      memcpy(out_ev_or_od+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out_ev_or_od);
  }
  THREADABLE_FUNCTION_END
  
  //paste the even and odd parts of a vector into a full lx vector
  THREADABLE_FUNCTION_3ARG(paste_eo_parts_into_lx_vector_internal, char*,out_lx, char**,in_eo, size_t,bps)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    //paste
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(eo,0,loc_volh)
	memcpy(out_lx+bps*loclx_of_loceo[par][eo],in_eo[par]+bps*eo,bps);
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out_lx);
  }
  THREADABLE_FUNCTION_END
  
  /////////////////////
  
  //remap using a certain local remapper
  THREADABLE_FUNCTION_4ARG(remap_lx_vector_internal, char*,out, char*,in, size_t,bps, int*,remap)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      memcpy(out+bps*remap[ivol],in+bps*ivol,bps);
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //remap eo using a certain remapper
  THREADABLE_FUNCTION_4ARG(remap_eo_vector_internal, char**,out, char**,in, size_t,bps, int**,remap)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivolh,0,loc_volh)
	memcpy(out[par]+bps*remap[par][ivolh],in[par]+bps*ivolh,bps);
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out[0]);
    set_borders_invalid(out[1]);
  }
  THREADABLE_FUNCTION_END
}
