#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <string.h>

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry_eo.hpp"
#include "geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  /////////////////////
  
  //remap using a certain local remapper
  void remap_vector_internal(char* out,char* in,size_t bps,int* dest_of_source,int length)
  {
    
    START_TIMING(remap_time,nremap);
    
    NISSA_PARALLEL_LOOP(i,0,length)
      memcpy(out+bps*dest_of_source[i],in+bps*i,bps);
    NISSA_PARALLEL_LOOP_END;
    
    STOP_TIMING(remap_time);
    
    set_borders_invalid(out);
  }
  
  ///////////////////////
  
  // //pack the 8 links attached to a given site
  // void remap_loceo_conf_to_Lebeo_oct(oct_su3* out,quad_su3** in,int par)
  // {
    
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
}
