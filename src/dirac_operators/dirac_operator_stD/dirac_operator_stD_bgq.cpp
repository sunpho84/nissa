#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "base/global_variables.h"
#include "base/thread_macros.h"
#include "base/vectors.h"
#include "bgq/bgq_macros.h"
#include "bgq/staggered_hopping_matrix_eo_or_oe_bgq.h"
#include "new_types/complex.h"
#ifdef USE_THREADS
 #include "routines/thread.h"
#endif

THREADABLE_FUNCTION_5ARG(apply_stD2ee_m2_bgq, bi_color*,bi_out, bi_oct_su3**,bi_conf, bi_color*,bi_temp, double,mass2, bi_color*,bi_in)
{
  //----------------------looping on E--------------------
  const int OE=0;
  
  //compute on the surface and start communications
  apply_staggered_hopping_matrix_oe_or_eo_bgq_nocomm_nobarrier(bi_conf,0,vsurf_volh,bi_in,OE);
  start_staggered_hopping_matrix_oe_or_eo_bgq_communications();
  
  //compute on the bulk and finish communications
  apply_staggered_hopping_matrix_oe_or_eo_bgq_nocomm_nobarrier(bi_conf,vsurf_volh,loc_volh/2,bi_in,OE);
  finish_staggered_hopping_matrix_oe_or_eo_bgq_communications(OE);

  //put the eight pieces together
  hopping_matrix_eo_or_eo_expand_to_D(bi_temp);
  
  //----------------------looping on O--------------------
  const int EO=1;  
  
  //compute on the surface and start communications
  apply_staggered_hopping_matrix_oe_or_eo_bgq_nocomm_nobarrier(bi_conf,0,vsurf_volh,bi_temp,EO);
  start_staggered_hopping_matrix_oe_or_eo_bgq_communications();
  
  //compute on the bulk and finish communications
  apply_staggered_hopping_matrix_oe_or_eo_bgq_nocomm_nobarrier(bi_conf,vsurf_volh,loc_volh/2,bi_temp,EO);
  finish_staggered_hopping_matrix_oe_or_eo_bgq_communications(EO);

  //put the eight pieces subtracting them from diag (in fact one of the two D is daggered)
  hopping_matrix_eo_or_eo_expand_to_D_subtract_from_mass2_times_in(bi_out,mass2,bi_in);
}}
