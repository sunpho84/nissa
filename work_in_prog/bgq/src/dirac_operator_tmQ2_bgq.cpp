#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

 #include <math.h>

#include "nissa.h"

#include "new_vars_and_types.h"
#include "bgq_macros.h"

/*
  In bgq version we merge to different virtual nodes along t directions,
  so that only site with time coordinate between 0 and T/2-1 must be considered.
  Since we want to load everything sequentially, we need to duplicate the gauge
  configuration.

  We apply first time the hopping matrix scanning on sink index, then store the 8 contributions
  to each source separately.
*/

THREADABLE_FUNCTION_3ARG(apply_bgq_Wilson_hopping_matrix, bi_halfspincolor*,out, bi_oct_su3*,conf, bi_spincolor*,in)
{
  GET_THREAD_ID();

  //check borders: they must all be filled
  if(IS_MASTER_THREAD)
    for(int ivol=0;ivol<loc_volh;ivol++)
      CHECK_BI_SPINCOLOR_DIFF_FROM_ZERO(in[ivol],ivol);
  
  NISSA_PARALLEL_LOOP(ibgqlx,0,loc_volh)
    {
      /*debug
	master_printf("ibgqlx=%d\n",ibgqlx);*/
      bi_spincolor temp;
      
      //take short access to link and output indexing
      int *ind=bgq_hopping_matrix_output_index+ibgqlx*8;
      bi_su3 *links=(bi_su3*)(conf+ibgqlx);
      
      //T backward derivative
      HOPMATR_TBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(out[ind[0]],links[0],temp);
      
      //X backward derivative
      HOPMATR_XBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(out[ind[1]],links[1],temp);
      
      //Y backward derivative
      HOPMATR_YBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(out[ind[2]],links[2],temp);
      
      //Z backward derivative
      HOPMATR_ZBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(out[ind[3]],links[3],temp);
      
      //T forward derivative
      HOPMATR_TFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR(out[ind[4]],links[4],temp);
      
      //X forward derivative
      HOPMATR_XFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR(out[ind[5]],links[5],temp);
      
      //Y forward derivative
      HOPMATR_YFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR(out[ind[6]],links[6],temp);
      
      //Z forward derivative
      HOPMATR_ZFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR(out[ind[7]],links[7],temp);
      
      /*debug
	for(int i=0;i<8;i++)
	master_printf("%d filling %04d\n",i,ind[i]);*/
    }
  
  /*debug 
    if(IS_MASTER_THREAD) //check borders: they must all be filled
    for(int ivol=0;ivol<bgqlx_vbord_vol;ivol++)
    CHECK_BI_SPINCOLOR_DIFF_FROM_ZERO(out[ivol],ivol);*/
  
  thread_barrier(HOPPING_MATRIX_APPLICATION_BARRIER);
  
  /////////////////////////////////// perform VN communications ////////////////////////////
  
  //split bw T border: 0 VN node goes to bw out border, 1 VN goes to VN 0
  
  
}}
