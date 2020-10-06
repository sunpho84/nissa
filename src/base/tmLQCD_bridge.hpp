#ifndef _TMLQCD_BRIDGE_HPP
#define _TMLQCD_BRIDGE_HPP

//include tmLQCD inside a clean namespace
namespace tmLQCD
{
  #include <tmLQCD.h>
}

#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"

#ifndef EXTERN_TMLQCD_BRIDGE
 #define EXTERN_TMLQCD_BRIDGE extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace nissa
{
  EXTERN_TMLQCD_BRIDGE quad_su3 *external_conf_to_tmLQCD_handle;
  EXTERN_TMLQCD_BRIDGE int use_tmLQCD INIT_TO(1);
  EXTERN_TMLQCD_BRIDGE int tmLQCD_initialized INIT_TO(0);
  
  //importing the finalizer
  using tmLQCD::tmLQCD_finalise;
  
  void tmLQCD_init();
  inline void nissa_feed_conf_to_tmLQCD(tmLQCD::su3 *ext,int t,int x,int y,int z,int mu)
  {su3_copy(*(su3*)(void*)ext,external_conf_to_tmLQCD_handle[loclx_of_coord_list(t,x,y,z)][mu]);}
  void export_gauge_conf_to_tmLQCD(quad_su3 *conf_lx);
  FILE* open_prepare_input_file_for_tmLQCD();
}

#undef EXTERN_TMLQCD_BRIDGE
#undef INIT_TO

#endif
