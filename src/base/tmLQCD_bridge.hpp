#ifndef _TMLQCD_BRIDGE_HPP
#define _TMLQCD_BRIDGE_HPP

#ifndef EXTERN_BRIDGE
 #define EXTERN_BRIDGE extern
#endif

//include tmLQCD inside a clean namespace
namespace tmLQCD
{
  #include <tmLQCD.h>
}

#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  EXTERN_BRIDGE quad_su3 *external_conf_to_tmLQCD_handle;
  
  //importing the finalizer
  using tmLQCD::tmLQCD_finalise;
  
  void tmLQCD_init();
  inline void nissa_feed_conf_to_tmLQCD(tmLQCD::su3 *ext,int t,int x,int y,int z,int mu)
  {su3_copy(*(su3*)(void*)ext,external_conf_to_tmLQCD_handle[loclx_of_coord_list(t,x,y,z)][mu]);}
  void export_gauge_conf_to_tmLQCD(quad_su3 *conf_lx);
  FILE* open_prepare_input_file_for_tmLQCD();
}

#undef EXTERN_BRIDGE

#endif
