#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry_eo.hpp"
#include "geometry_lx.hpp"
#include "geometry_Leb.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //separate the even and odd part of a vector
  THREADABLE_FUNCTION_3ARG(split_lx_vector_into_eo_parts_internal, eo_ptr<void>,out_eo, void*,in_lx, size_t,bps)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    //split
    NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
      memcpy((char*)(out_eo[loclx_parity[loclx]])+bps*loceo_of_loclx[loclx],(char*)in_lx+bps*loclx,bps);
    NISSA_PARALLEL_LOOP_END;
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out_eo[0]);
    set_borders_invalid(out_eo[1]);
  }
  THREADABLE_FUNCTION_END
  
  //separate the even and odd part of a vector
  THREADABLE_FUNCTION_4ARG(get_evn_or_odd_part_of_lx_vector_internal, void*,out_ev_or_od, void*,in_lx, size_t,bps, int,par)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    //get
    NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
      if(loclx_parity[loclx]==par)
	memcpy((char*)out_ev_or_od+bps*loceo_of_loclx[loclx],(char*)in_lx+bps*loclx,bps);
    NISSA_PARALLEL_LOOP_END;
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out_ev_or_od);
  }
  THREADABLE_FUNCTION_END
  
  //paste the even and odd parts of a vector into a full lx vector
  THREADABLE_FUNCTION_3ARG(paste_eo_parts_into_lx_vector_internal, void*,out_lx, eo_ptr<void>,in_eo, size_t,bps)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    //paste
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(eo,0,loc_volh)
	memcpy((char*)out_lx+bps*loclx_of_loceo[par][eo],(char*)(in_eo[par])+bps*eo,bps);
    NISSA_PARALLEL_LOOP_END;
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out_lx);
  }
  THREADABLE_FUNCTION_END
  
  /////////////////////
  
  //remap using a certain local remapper
  THREADABLE_FUNCTION_5ARG(remap_vector_internal, char*,out, char*,in, size_t,bps, int*,dest_of_source, int,length)
  {
    GET_THREAD_ID();
    
    START_TIMING(remap_time,nremap);
    
    NISSA_PARALLEL_LOOP(i,0,length)
      memcpy(out+bps*dest_of_source[i],in+bps*i,bps);
    NISSA_PARALLEL_LOOP_END;
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  ///////////////////////
  
  // //pack the 8 links attached to a given site
  // THREADABLE_FUNCTION_3ARG(remap_loceo_conf_to_Lebeo_oct, oct_su3*,out, quad_su3**,in, int,par)
  // {
  //   GET_THREAD_ID();
    
  //   START_TIMING(remap_time,nremap);
    
  //   communicate_ev_and_od_quad_su3_borders(in);
    
  //   NISSA_PARALLEL_LOOP(Lebdest,0,loc_volh)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	{
  // 	  su3_copy(out[Lebdest][NDIM+mu],in[par][loceo_of_Lebeo[par][Lebdest]][mu]);
  // 	  su3_copy(out[Lebdest][mu],in[!par][loceo_of_Lebeo[!par][Lebeo_neighdw[par][Lebdest][mu]]][mu]);
  // 	}
  //   NISSA_PARALLEL_LOOP_END;
    
  //   set_borders_invalid(out);
    
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
}
