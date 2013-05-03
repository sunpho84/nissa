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
  
  //reset temp with -1, nissa_recv_buf with -2 and nissa_send_buf with -3
  NISSA_PARALLEL_LOOP(ivol,0,loc_volh*8)
    ((bi_halfspincolor*)temp)[ivol][0][0][0][0]=-1;
  NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
    {
      ((bi_halfspincolor*)nissa_recv_buf)[ibord][0][0][0][0]=-2;
      ((bi_halfspincolor*)nissa_send_buf)[ibord][0][0][0][0]=-3;
    }
  
  ////////////////////////////////////////// core /////////////////////////////////////////
  
  apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(conf,0,bgq_vsurf_vol,in);
  start_bgq_Wilson_hopping_matrix_communications();
  apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(conf,bgq_vsurf_vol,loc_volh,in);
  finish_bgq_Wilson_hopping_matrix_communications();

  ////////////////////////// these operations are done only temporary /////////////////////
  
  int io=0;
  int ivn=io>=loc_volh;
  int ib=bgqlx_of_loclx[io%loc_volh];
  master_printf("%d %d %d:\n",io,ivn,ib);
  for(int i=0;i<8;i++)
    master_printf("i %d, %16.16lg\n",i,temp[ib*8+i][0][0][ivn][0]);
  
  //summ the 8 components of the output
  thread_barrier(HOPPING_MATRIX_APPLICATION_BARRIER);
  NISSA_PARALLEL_LOOP(X,0,loc_volh)
    {
      bi_spincolor t;
      memset(t,0,sizeof(bi_spincolor));
      for(int i=0;i<8;i++)
	BI_HALFSPINCOLOR_SUMMASSIGN(t,temp[X*8+i]);
      BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(out[loclx_of_bgqlx[X]],out[loclx_of_bgqlx[X]+loc_volh],t);
    }
  
  thread_barrier(HOPPING_MATRIX_APPLICATION_BARRIER);
  
  master_printf("tot: %lg\n",out[io][0][0][0]);
  
  //unbind the pointer
  nissa_free(bgq_hopping_matrix_output_pointer);
  nissa_free(bgq_hopping_matrix_output_T_buffer);
  
  //free temp buffer
  nissa_free(temp);
}}
