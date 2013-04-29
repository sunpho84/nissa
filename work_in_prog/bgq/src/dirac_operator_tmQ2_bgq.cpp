#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

 #include <math.h>

#include "nissa.h"

#include "bgq_macros.h"
#include "geometry_bgq.h"
#include "hopping_matrix_bgq.h"
#include "new_vars_and_types.h"

/*
  application of hopping matrix goes among the following steps:
   1) hopping matrix on the surf
   2) start communications
   3) hopping matrix on the bulk
   4) finish communications
*/
THREADABLE_FUNCTION_5ARG(apply_tmQ_bgq, spincolor*,out, bi_oct_su3*,conf, double,kappa, double,mu, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  //allocate a temporary vector to apply hopping matrix
  bi_halfspincolor *temp=nissa_malloc("temp",8*loc_volh,bi_halfspincolor);
  
  //define output index pointer binded to temp
  define_bgq_hopping_matrix_lx_output_pointers_and_T_buffers(temp);
  
  ////////////////////////////////////////// core /////////////////////////////////////////
  
  apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(conf,0,bgq_vsurf_vol,in);
  start_bgq_Wilson_hopping_matrix_communications();
  apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(conf,bgq_vsurf_vol,loc_volh,in);
  finish_bgq_Wilson_hopping_matrix_communications();
  
  ////////////////////////// these operations are done only temporary /////////////////////
  
  //summ the 8 components of the output
  NISSA_PARALLEL_LOOP(X,0,loc_volh)
    {
      bi_spincolor t;
      memset(t,0,sizeof(bi_spincolor));
      for(int i=0;i<8;i++)
	BI_HALFSPINCOLOR_SUMMASSIGN(t,temp[X*8+i]);
      BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(out[loclx_of_bgqlx[X]],out[loclx_of_bgqlx[X]+loc_volh],t);
    }
  
  //unbind the pointer
  nissa_free(bgq_hopping_matrix_output_pointer);
  nissa_free(bgq_hopping_matrix_output_T_buffer);
  
  //free temp buffer
  nissa_free(temp);
}}
