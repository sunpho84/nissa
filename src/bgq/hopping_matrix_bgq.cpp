#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../base/global_variables.h"
#include "../base/thread_macros.h"
#include "../communicate/borders.h"
#include "../new_types/complex.h"
#include "../new_types/new_types_definitions.h"
#ifdef USE_THREADS
 #include "../routines/thread.h"
#endif

#include "bgq_macros.h"

/*
  In bgq version we merge two sites along t directions, that is, (t,x,y,z) and (t+T/2,x,y,z),
  so that only site with time coordinate between 0 and T/2-1 must be considered.
  Since we want to load everything sequentially, we need to duplicate the gauge configuration.

  We apply hopping matrix scanning on sink index, then store the 8 contributions to each source separately.
  They can be summed outside according to external usage.
  
  We do not sync and do not perform any communication, so this must be done outside.
*/

#define PROJ_HEADER(A)				\
  REORDER_BARRIER();				\
  CACHE_PREFETCH(out+A);			\
  BI_SU3_PREFETCH_NEXT(links[A])

#ifdef BGQ
 #define BORDER_SITE_COPY(base_out,base_in,isrc)	\
  {							\
    void *in=base_in+isrc;				\
    CACHE_PREFETCH(base_out+isrc);			\
    DECLARE_REG_BI_HALFSPINCOLOR(reg_temp);		\
    BI_HALFSPINCOLOR_PREFETCH_NEXT(in);			\
    REG_LOAD_BI_HALFSPINCOLOR(reg_temp,in);		\
    void *out=base_out[isrc];				\
    STORE_REG_BI_HALFSPINCOLOR(out,reg_temp);		\
  }
#else
 #define BORDER_SITE_COPY(base_out,base_in,isrc) BI_HALFSPINCOLOR_COPY((*(base_out[isrc])),base_in[isrc]);
#endif

THREADABLE_FUNCTION_4ARG(apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier, bi_oct_su3*,conf, int,istart, int,iend, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  NISSA_PARALLEL_LOOP(ibgqlx,istart,iend)
    {
      //take short access to link and output indexing
      bi_halfspincolor **out=bgq_hopping_matrix_output_pointer+ibgqlx*8;
      bi_su3 *links=(bi_su3*)(conf+ibgqlx);
      
      //declare
      DECLARE_REG_BI_SPINCOLOR(reg_in);
      DECLARE_REG_BI_HALFSPINCOLOR(reg_proj);
      
      //load in
      REG_LOAD_BI_SPINCOLOR(reg_in,in[ibgqlx]);
      
      //T backward scatter (forward derivative)
      PROJ_HEADER(0);
      REG_BI_COLOR_SUMM(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_SUMM(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[0])),links[0],reg_proj);
      
      //X backward scatter (forward derivative)
      PROJ_HEADER(1);
      REG_BI_COLOR_ISUMM(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_ISUMM(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[1])),links[1],reg_proj);
      
      //Y backward scatter (forward derivative)
      PROJ_HEADER(2);
      REG_BI_COLOR_SUMM(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_SUBT(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[2])),links[2],reg_proj);
      
      //Z backward scatter (forward derivative)
      PROJ_HEADER(3);
      REG_BI_COLOR_ISUMM(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_ISUBT(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[3])),links[3],reg_proj);
      
      //T forward scatter (backward derivative)
      PROJ_HEADER(4);
      REG_BI_COLOR_SUBT(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_SUBT(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[4])),links[4],reg_proj);
      
      //X forward scatter (backward derivative)
      PROJ_HEADER(5);
      REG_BI_COLOR_ISUBT(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_ISUBT(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[5])),links[5],reg_proj);

      //Y forward scatter (backward derivative)
      PROJ_HEADER(6);
      REG_BI_COLOR_SUBT(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_SUMM(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[6])),links[6],reg_proj);
      
      //Z forward scatter (backward derivative)
      PROJ_HEADER(7);
      REG_BI_COLOR_ISUBT(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_ISUMM(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE((*(out[7])),links[7],reg_proj);
    }
}}
  
//swap border between VN and, if T is really parallelized, fill send buffers
THREADABLE_FUNCTION_0ARG(bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill)
{
  GET_THREAD_ID();
  
  //T buffers might have been filled in another order
  THREAD_BARRIER();
  
  //form two thread team
  FORM_TWO_THREAD_TEAMS();
  
  ///////////////////////// bw scattered T derivative (fw derivative)  ////////////////////////
  
  if(is_in_first_team)
    if(paral_dir[0])
      {
	//split bw T border: VN 0 goes to bw out border (first half), VN 1 goes to VN 0 fw surf
	NISSA_CHUNK_LOOP(isrc,0,bgqlx_t_vbord_vol/2,thread_in_team_id,nthreads_in_team)
	  {
	    //non-local shuffling: must enter bw buffer for direction 0
	    int idst_buf=isrc;
#ifndef EXP_BGQ
	    //local shuffling: inside the result border buffer, data is stored as in the border, so isrc access to t-face
	    //summing loc_volh we get to the loc_size[0]/2 plane, and its backward neighbour is the destination
	    int idst_loc=8*bgqlx_of_loclx[loclx_neighdw[isrc+loc_volh][0]];
	    halfspincolor temp;
	    BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(((halfspincolor*)nissa_send_buf)[idst_buf],       //take VN=0
					      temp,bgq_hopping_matrix_output_T_buffer[isrc]);
	    HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_data[idst_loc],temp,0); //copy VN=1 to VN=0
#else
	    //load the first
	    DECLARE_REG_BI_HALFSPINCOLOR(in0);
	    REG_LOAD_BI_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_T_buffer[isrc]);
	    //load the second
	    DECLARE_REG_BI_HALFSPINCOLOR(in1);
	    REG_LOAD_BI_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_T_buffer[isrc+1]);
	    //merge the two and save
	    DECLARE_REG_BI_HALFSPINCOLOR(to_buf);
	    REG_BI_HALFSPINCOLOR_V0_MERGE(to_buf,in0,in1);
	    STORE_REG_BI_HALFSPINCOLOR(((halfspincolor*)nissa_send_buf)[idst_buf],to_buf);
	    
	    //paired
	    isrc++;
#endif
	  }
      }
    else
      //we have only to transpose between VN, and no real communication happens
      NISSA_CHUNK_LOOP(isrc,0,bgqlx_t_vbord_vol/2,thread_in_team_id,nthreads_in_team)
	{
	  //look at previous idst_loc for explenation
	  int idst=8*bgqlx_of_loclx[loclx_neighdw[isrc+loc_volh][0]];
	  BI_HALFSPINCOLOR_TRANSPOSE(bgq_hopping_matrix_output_data[idst],bgq_hopping_matrix_output_T_buffer[isrc]);
	}
  
  ///////////////////////// fw scattered T derivative (bw derivative)  ////////////////////////
  
  if(is_in_second_team)
    if(paral_dir[0])
      {
	//split fw T border: VN 1 goes to fw out border (second half), VN 0 goes to VN 1 of t=0 surf
	NISSA_CHUNK_LOOP(base_isrc,0,bgqlx_t_vbord_vol/2,thread_in_team_id,nthreads_in_team)
	  {
	    //the source starts at the middle of result border buffer
	    int isrc=base_isrc+bgqlx_t_vbord_vol/2;
	    //non-local shuffling: must enter fw buffer for direction 0
	    int idst_buf=bord_volh+base_isrc;
#ifndef EXP_BGQ
	    //local shuffling: inside the result border buffer, data is stored as in the border, so base_isrc
	    //points to the correct destination, but we must shift for terms of bw der
	    int idst_loc=8*base_isrc+4;
	    halfspincolor temp;
	    BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(temp,((halfspincolor*)nissa_send_buf)[idst_buf],
					      bgq_hopping_matrix_output_T_buffer[isrc]);
	    HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_data[idst_loc],temp,1); //copy VN=0 to VN=1
#else
	    //load the first
	    DECLARE_REG_BI_HALFSPINCOLOR(in0);
	    REG_LOAD_BI_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_T_buffer[isrc]);
	    //load the second
	    DECLARE_REG_BI_HALFSPINCOLOR(in1);
	    REG_LOAD_BI_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_T_buffer[isrc+1]);
	    //merge the two and save
	    DECLARE_REG_BI_HALFSPINCOLOR(to_buf);
	    REG_BI_HALFSPINCOLOR_V1_MERGE(to_buf,in0,in1);
	    STORE_REG_BI_HALFSPINCOLOR(((halfspincolor*)nissa_send_buf)[idst_buf],to_buf);
	    
	    //paired
	    base_isrc++;
#endif
	  }
      }
    else
      //we have only to transpose between VN
      NISSA_CHUNK_LOOP(base_isrc,0,bgqlx_t_vbord_vol/2,thread_in_team_id,nthreads_in_team)
	{
	  //the source starts at the middle of result border buffer
	  int isrc=base_isrc+bgqlx_t_vbord_vol/2;
	  //look at previous VN=0 for explenation
	  int idst=8*(bord_volh+base_isrc)+4;
	  BI_HALFSPINCOLOR_TRANSPOSE(bgq_hopping_matrix_output_data[idst],bgq_hopping_matrix_output_T_buffer[isrc]);
	}

  
  ///////////// debugging 
#if 0
  thread_barrier();
  for(int ivol=0;ivol<bord_vol;ivol++)
    for(int idouble=0;idouble<sizeof(halfspincolor)/sizeof(double);idouble++)
      {
	double a=((double*)(((halfspincolor*)nissa_send_buf)[ivol]))[idouble];
	if(isnan(a)) crash("idouble %d ivol %d",idouble,ivol);
      }
#endif
}}

//perform communications between VN and start all the communications between nodes
THREADABLE_FUNCTION_0ARG(start_Wilson_hopping_matrix_bgq_binded_communications)
{
  //shuffle data between virtual nodes and fill T out buffer
  bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill();
  
  //after the barrier, all buffers are filled and communications can start
  THREAD_BARRIER();
  
  //start communications of scattered data to other nodes
  comm_start(lx_halfspincolor_comm);
}}

//finish the communications and put in place the communicated data
THREADABLE_FUNCTION_0ARG(finish_Wilson_hopping_matrix_bgq_binded_communications)
{
  GET_THREAD_ID();
  FORM_TWO_THREAD_TEAMS();
  
  //wait communications end
  comm_wait(lx_halfspincolor_comm);
  
  //T bw border (bw derivative): data goes to VN 0
  if(is_in_first_team)
    {
      bi_halfspincolor **base_out=bgq_hopping_matrix_final_output;
      halfspincolor *base_in=(halfspincolor*)nissa_recv_buf;
      
#ifndef EXP_BGQ
      NISSA_CHUNK_LOOP(isrc,0,bord_dir_vol[0],thread_in_team_id,nthreads_in_team)
	HALFSPINCOLOR_TO_BI_HALFSPINCOLOR((*(base_out[isrc])),base_in[isrc],0);
#else
      NISSA_CHUNK_LOOP(isrc,0,bord_dir_vol[0],thread_in_team_id,nthreads_in_team)
	{
	  //VN=0 must be filled with border
	  DECLARE_REG_BI_HALFSPINCOLOR(in0);
	  REG_LOAD_BI_HALFSPINCOLOR(in0,base_in[isrc]);
	  //VN=1 with buf0
	  DECLARE_REG_BI_HALFSPINCOLOR(in1);
	  REG_LOAD_BI_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_T_buffer[isrc+bgqlx_t_vbord_vol/2]);
	  //merge and save
	  DECLARE_REG_BI_HALFSPINCOLOR(to_dest);
	  REG_BI_HALFSPINCOLOR_V0_MERGE(to_dest,in0,in1);
	  STORE_REG_BI_HALFSPINCOLOR((*(base_out[isrc])),to_dest);
	  
	  //paired
	  isrc++;
	  
	  //VN=1 with buf1
	  REG_LOAD_BI_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_T_buffer[isrc+bgqlx_t_vbord_vol/2]);
	  //merge and save
	  REG_BI_HALFSPINCOLOR_V10_MERGE(to_dest,in0,in1);
	  STORE_REG_BI_HALFSPINCOLOR((*(base_out[isrc])),to_dest);
	}
#endif
    }
  
  //other 3 bw borders
  if(is_in_first_team)
    {
      bi_halfspincolor **base_out=bgq_hopping_matrix_final_output+bord_dir_vol[0];
      bi_halfspincolor *base_in=(bi_halfspincolor*)(nissa_recv_buf+sizeof(bi_halfspincolor)*bord_dir_vol[0]/2);
      NISSA_CHUNK_LOOP(isrc,0,bord_vol/4-bord_dir_vol[0]/2,thread_in_team_id,nthreads_in_team)
	BORDER_SITE_COPY(base_out,base_in,isrc);
    }
  
  //T fw border (fw derivative): data goes to VN 1
  if(is_in_second_team)
    {
      bi_halfspincolor **base_out=bgq_hopping_matrix_final_output+bord_dir_vol[0]/2+bord_vol/4;
      halfspincolor *base_in=(halfspincolor*)(nissa_recv_buf+sizeof(bi_halfspincolor)*bord_vol/4);
#ifndef EXP_BGQ
      NISSA_CHUNK_LOOP(isrc,0,bord_dir_vol[0],thread_in_team_id,nthreads_in_team)
	HALFSPINCOLOR_TO_BI_HALFSPINCOLOR((*(base_out[isrc])),base_in[isrc],1);
#else
      NISSA_CHUNK_LOOP(isrc,0,bord_dir_vol[0],thread_in_team_id,nthreads_in_team)
	{
	  //VN=0 with buf1
	  DECLARE_REG_BI_HALFSPINCOLOR(in0);
	  REG_LOAD_BI_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_T_buffer[isrc]);
	  //VN=1 must be filled with border 0
	  DECLARE_REG_BI_HALFSPINCOLOR(in1);
	  REG_LOAD_BI_HALFSPINCOLOR(in1,base_in[isrc]);
	  //merge and save
	  DECLARE_REG_BI_HALFSPINCOLOR(to_dest);
	  REG_BI_HALFSPINCOLOR_V10_MERGE(to_dest,in0,in1);
	  STORE_REG_BI_HALFSPINCOLOR((*(base_out[isrc])),to_dest);
	  
	  //paired
	  isrc++;
	  
	  //VN=0 with buf1
	  REG_LOAD_BI_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_T_buffer[isrc]);
	  //merge and save
	  REG_BI_HALFSPINCOLOR_V1_MERGE(to_dest,in0,in1);
	  STORE_REG_BI_HALFSPINCOLOR((*(base_out[isrc])),to_dest);
	}
#endif
    }
  
  //other 3 fw borders
  if(is_in_second_team)
    {
      bi_halfspincolor **base_out=bgq_hopping_matrix_final_output+3*bord_dir_vol[0]/2+bord_vol/4;
      bi_halfspincolor *base_in=(bi_halfspincolor*)(nissa_recv_buf+sizeof(bi_halfspincolor)*(bord_dir_vol[0]/2+bord_vol/4));
      NISSA_CHUNK_LOOP(isrc,0,bord_vol/4-bord_dir_vol[0]/2,thread_in_team_id,nthreads_in_team)
	BORDER_SITE_COPY(base_out,base_in,isrc);
    }
}}
