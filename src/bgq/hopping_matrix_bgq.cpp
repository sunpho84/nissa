#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../base/global_variables.h"
#include "../base/thread_macros.h"
#include "../communicate/borders.h"
#include "../new_types/complex.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/thread.h"

#include "bgq_macros.h"

/*
  In bgq version we merge two sites along t directions, that is, (t,x,y,z) and (t+T/2,x,y,z),
  so that only site with time coordinate between 0 and T/2-1 must be considered.
  Since we want to load everything sequentially, we need to duplicate the gauge configuration.

  We apply hopping matrix scanning on sink index, then store the 8 contributions to each source separately.
  They can be summed outside according to external usage.
  
  We do not sync and do not perform any communication, so this must be done outside.
*/

THREADABLE_FUNCTION_4ARG(apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier, bi_oct_su3*,conf, int,istart, int,iend, bi_spincolor*,in)
{
  GET_THREAD_ID();

  NISSA_PARALLEL_LOOP(ibgqlx,istart,iend)
    {
      //take short access to link and output indexing
      bi_halfspincolor **out=bgq_hopping_matrix_output_pointer+ibgqlx*8;
      bi_su3 *links=(bi_su3*)(conf+ibgqlx);
      
      bi_halfspincolor temp;
      
      //T backward scatter (forward derivative)
      HOPMATR_TFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR((*(out[0])),links[0],temp);
      
      //X backward scatter (forward derivative)
      HOPMATR_XFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR((*(out[1])),links[1],temp);
      
      //Y backward scatter (forward derivative)
      HOPMATR_YFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR((*(out[2])),links[2],temp);
      
      //Z backward scatter (forward derivative)
      HOPMATR_ZFW_PROJ(temp,in[ibgqlx]);
      BI_SU3_PROD_BI_HALFSPINCOLOR((*(out[3])),links[3],temp);
      
      //T forward scatter (backward derivative)
      HOPMATR_TBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR((*(out[4])),links[4],temp);
      
      //X forward scatter (backward derivative)
      HOPMATR_XBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR((*(out[5])),links[5],temp);
      
      //Y forward scatter (backward derivative)
      HOPMATR_YBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR((*(out[6])),links[6],temp);
      
      //Z forward scatter (backward derivative)
      HOPMATR_ZBW_PROJ(temp,in[ibgqlx]);
      BI_SU3_DAG_PROD_BI_HALFSPINCOLOR((*(out[7])),links[7],temp);
    }
}}
  
//swap border between VN and, if T is really parallelized, fill send buffers
void bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill()
{
  GET_THREAD_ID();
  
  //T buffers might have been filled in another order
  THREAD_BARRIER();
  
  ///////////////////////// bw scattered T derivative (fw derivative)  ////////////////////////
  
  if(paral_dir[0])
    {
      //split bw T border: VN 0 goes to bw out border (first half), VN 1 goes to VN 0 fw surf
      NISSA_PARALLEL_LOOP(isrc,0,bgqlx_t_vbord_vol/2)
	{
	  //local shuffling: inside the result border buffer, data is stored as in the border, so isrc access to t- face
	  //summing loc_volh we get to the loc_size[0]/2 plane, and its backward neighbour is the destination
	  int idst_loc=8*bgqlx_of_loclx[loclx_neighdw[isrc+loc_volh][0]];
	  //non-local shuffling: must enter bw buffer for direction 0
	  int idst_buf=isrc;
	  halfspincolor temp;
	  BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(((halfspincolor*)nissa_send_buf)[idst_buf],        //take VN=0
	  				    temp,bgq_hopping_matrix_output_T_buffer[isrc]);
	  HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_data[idst_loc],temp,0); //copy VN=1 to VN=0
	}
    }
  else
    //we have only to transpose between VN, and no real communication happens
    NISSA_PARALLEL_LOOP(isrc,0,bgqlx_t_vbord_vol/2)
      {
	//look at previous idst_loc for explenation
	int idst=8*bgqlx_of_loclx[loclx_neighdw[isrc+loc_volh][0]];
	BI_HALFSPINCOLOR_TRANSPOSE(bgq_hopping_matrix_output_data[idst],bgq_hopping_matrix_output_T_buffer[isrc]);
      }
  
  ///////////////////////// fw scattered T derivative (bw derivative)  ////////////////////////
  
  if(paral_dir[0])
    {
      //split fw T border: VN 1 goes to fw out border (second half), VN 0 goes to VN 1 of t=0 surf
      NISSA_PARALLEL_LOOP(base_isrc,0,bgqlx_t_vbord_vol/2)
	{
	  //the source starts at the middle of result border buffer
	  int isrc=base_isrc+bgqlx_t_vbord_vol/2;
	  //local shuffling: inside the result border buffer, data is stored as in the border, so base_isrc
	  //points to the correct destination, but we must shift for terms of bw der
	  int idst_loc=8*base_isrc+4;
	  //non-local shuffling: must enter fw buffer for direction 0
	  int idst_buf=bord_volh+base_isrc;
	  halfspincolor temp;
	  
	  BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(temp,((halfspincolor*)nissa_send_buf)[idst_buf],
					    bgq_hopping_matrix_output_T_buffer[isrc]);
	  HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_data[idst_loc],temp,1); //copy VN=0 to VN=1
	}
    }
  else
    //we have only to transpose between VN
    NISSA_PARALLEL_LOOP(base_isrc,0,bgqlx_t_vbord_vol/2)
      {
	//the source starts at the middle of result border buffer
	int isrc=base_isrc+bgqlx_t_vbord_vol/2;
	//look at previous VN=0 for explenation
	int idst=8*(bord_volh+base_isrc)+4;
	BI_HALFSPINCOLOR_TRANSPOSE(bgq_hopping_matrix_output_data[idst],bgq_hopping_matrix_output_T_buffer[isrc]);
      }
}

//perform communications between VN and start all the communications between nodes
void start_Wilson_hopping_matrix_bgq_binded_communications()
{
  //shuffle data between virtual nodes and fill T out buffer
  bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill();
  
  //after the barrier, all buffers are filled and communications can start
  THREAD_BARRIER();
  
  //start communications of scattered data to other nodes
  comm_start(lx_halfspincolor_comm);
}

//finish the communications and put in place the communicated data
void finish_Wilson_hopping_matrix_bgq_binded_communications()
{
  GET_THREAD_ID();
  
  //wait communications end
  comm_wait(lx_halfspincolor_comm);
  
  //T bw border (bw derivative): data goes to VN 0
  if(paral_dir[0])
    NISSA_PARALLEL_LOOP(isrc,0,bgqlx_t_vbord_vol/2)
      {
	int idst=8*bgqlx_of_loclx[isrc]+4; //bw derivative contributions start from 4
	HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_data[idst],((halfspincolor*)nissa_recv_buf)[isrc],0);
      }
  
  //other 3 bw borders
  for(int mu=1;mu<4;mu++)
    if(paral_dir[mu])
      NISSA_PARALLEL_LOOP(base_isrc,0,bord_dir_vol[mu]/2)
	{
	  int isrc=bord_offset[mu]/2+base_isrc;
	  int idst=8*bgqlx_of_loclx[surflx_of_bordlx[bord_offset[mu]+base_isrc]]+4+mu;
	  BI_HALFSPINCOLOR_COPY(bgq_hopping_matrix_output_data[idst],((bi_halfspincolor*)nissa_recv_buf)[isrc]);
      }
  
  //T fw border (fw derivative): data goes to VN 1
  if(paral_dir[0])
    NISSA_PARALLEL_LOOP(base_isrc,0,bgqlx_t_vbord_vol/2)
      {
	int isrc=bord_volh+base_isrc; //we are considering bi_halfspincolor as 2 halfspincolor, so full bord_vol
	int idst=8*bgqlx_of_loclx[loclx_neighdw[base_isrc+loc_volh][0]]; //fw derivative contributions start from 0
	HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_data[idst],((halfspincolor*)nissa_recv_buf)[isrc],1);
      }
  
  //other 3 fw borders
  for(int mu=1;mu<4;mu++)
    if(paral_dir[mu])
      NISSA_PARALLEL_LOOP(base_isrc,0,bord_dir_vol[mu]/2)
	{
	  int isrc=bord_volh/2+bord_offset[mu]/2+base_isrc;
	  int idst=8*bgqlx_of_loclx[surflx_of_bordlx[bord_volh+bord_offset[mu]+base_isrc]]+mu;
	  BI_HALFSPINCOLOR_COPY(bgq_hopping_matrix_output_data[idst],((bi_halfspincolor*)nissa_recv_buf)[isrc]);
      }
}
