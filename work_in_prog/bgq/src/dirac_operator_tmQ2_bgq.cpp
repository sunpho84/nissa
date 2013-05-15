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
  then the application of Q requires to expand the halfspincolor into full
  spincolor and summ the diagonal part
*/

THREADABLE_FUNCTION_4ARG(hopping_matrix_expand_to_Q_and_summ_diag_term_bgq_binded, bi_spincolor*,out, double,kappa, double,mu, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  double A=-1/kappa,B=-2*mu;
  bi_complex diag[2]={{{+A,B},{+A,B}},{{-A,B},{-A,B}}};
  
  //wait that all the terms are put in place
  THREAD_BARRIER();
  
  NISSA_PARALLEL_LOOP(i,0,loc_volh)
    {
      bi_halfspincolor *piece=bgq_hopping_matrix_output_binded+i*8;
      bi_spincolor temp;
      
      //multiply in by the diagonal part
      DIAG_TMQ(temp,diag,in[i]);
      
      //summ 8 contributions
      TFW_DER_TMQ_EXP(temp,piece[0]);
      XFW_DER_TMQ_EXP(temp,piece[1]);
      YFW_DER_TMQ_EXP(temp,piece[2]);
      ZFW_DER_TMQ_EXP(temp,piece[3]);
      TBW_DER_TMQ_EXP(temp,piece[4]);
      XBW_DER_TMQ_EXP(temp,piece[5]);
      YBW_DER_TMQ_EXP(temp,piece[6]);
      ZBW_DER_TMQ_EXP(temp,piece[7]);

      //put final 0.5
      BI_SPINCOLOR_PROD_DOUBLE(out[i],temp,-0.5);
    }
}}

THREADABLE_FUNCTION_5ARG(apply_tmQ_bgq, bi_spincolor*,out, bi_oct_su3*,conf, double,kappa, double,mu, bi_spincolor*,in)
{
  //allocate a temporary vector to apply hopping matrix
  bi_halfspincolor *temp=nissa_malloc("temp",8*loc_volh,bi_halfspincolor);
  
  //define output index pointer binded to temp
  define_bgq_hopping_matrix_lx_output_pointers_and_T_buffers(temp);
  
  /*debug
  GET_THREAD_ID();
  //reset temp with -1, nissa_recv_buf with -2 and nissa_send_buf with -3
  NISSA_PARALLEL_LOOP(ivol,0,loc_volh*8)
  ((bi_halfspincolor*)temp)[ivol][0][0][0][0]=-1;
  NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
  {
  ((bi_halfspincolor*)nissa_recv_buf)[ibord][0][0][0][0]=-2;
  ((bi_halfspincolor*)nissa_send_buf)[ibord][0][0][0][0]=-3;
  }
  THREAD_BARRIER();
  */  
  ////////////////////////////////////////// core /////////////////////////////////////////
  
  apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(conf,0,bgq_vsurf_vol,in);
  start_Wilson_hopping_matrix_bgq_binded_communications();
  apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(conf,bgq_vsurf_vol,loc_volh,in);
  finish_Wilson_hopping_matrix_bgq_binded_communications();
  hopping_matrix_expand_to_Q_and_summ_diag_term_bgq_binded(out,kappa,mu,in);
  
  ////////////////////////// these operations are done only temporary /////////////////////
  
  /*non working debug
    int io=0;
    int ivn=io>=loc_volh;
    int ib=bgqlx_of_loclx[io%loc_volh];
    master_printf("%d %d %d:\n",io,ivn,ib);
    for(int i=0;i<8;i++)
    master_printf("i %d, %16.16lg\n",i,temp[ib*8+i][0][0][ivn][0]);
    
    //summ the 8 components of the output
    THREAD_BARRIER();
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
    {
    bi_spincolor t;
    memset(t,0,sizeof(bi_spincolor));
    for(int i=0;i<8;i++)
    BI_HALFSPINCOLOR_SUMMASSIGN(t,temp[X*8+i]);
    BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(out[loclx_of_bgqlx[X]],out[loclx_of_bgqlx[X]+loc_volh],t);
    }
    
    THREAD_BARRIER();
    
    master_printf("tot: %lg\n",out[io][0][0][0]);
  */

  //unbind the pointer
  nissa_free(bgq_hopping_matrix_output_pointer);
  nissa_free(bgq_hopping_matrix_output_T_buffer);
  
  //free temp buffer
  nissa_free(temp);
}}
