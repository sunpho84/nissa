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

THREADABLE_FUNCTION_4ARG(apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier, bi_oct_su3*,conf, int,istart, int,iend, bi_spincolor*,in)
{
  GET_THREAD_ID();

  NISSA_PARALLEL_LOOP(ibgqlx,istart,iend)
    {
      //take short access to link and output indexing
      bi_halfspincolor **out=bgq_hopping_matrix_output_pointer+ibgqlx*8;
      bi_su3 *links=(bi_su3*)(conf+ibgqlx);
      
      bi_spincolor temp;
      
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
  
//swap border between VN and if T is paral fill send buffers
void bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill()
{
  GET_THREAD_ID();
  
  //T buffers might be filled in another order
  thread_barrier(HOPPING_MATRIX_APPLICATION_BARRIER);
  
  ///////////////////////// bw scattered T derivative (fw derivative)  ////////////////////////
  
  if(paral_dir[0])
    {
      //split bw T border: VN 0 goes to bw out border, VN 1 goes to VN 0
      NISSA_PARALLEL_LOOP(isrc,0,bgqlx_t_vbord_vol/2)
	{
	  //local shuffling: inside the result border buffer, data is stored as in the border, so isrc access to t- face
	  //summing loc_volh we get to the loc_size[0]/2 plane, and its backward neighbour is the destination
	  int idst_loc=8*bgqlx_of_loclx[loclx_neighdw[isrc+loc_volh][0]];
	  //non-local shuffling: must enter bw buffer for direction 0
	  int idst_buf=bord_volh+isrc;
	  halfspincolor temp;
	  BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(((halfspincolor*)nissa_send_buf)[idst_buf], //take VN=0
					    temp,bgq_hopping_matrix_output_T_buffer[isrc]);
	  HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_binded[idst_loc],temp,0); //copy VN=1 to VN=0
	}
    }
  else
    //we have only to transpose between VN, and no real communication happens
    NISSA_PARALLEL_LOOP(isrc,0,bgqlx_t_vbord_vol/2)
      {
	//look at previous idst_loc for explanation
	int idst=8*bgqlx_of_loclx[loclx_neighdw[isrc+loc_volh][0]];
	BI_HALFSPINCOLOR_TRANSPOSE(bgq_hopping_matrix_output_binded[idst],bgq_hopping_matrix_output_T_buffer[isrc]);
      }
  
  ///////////////////////// fw scattered T derivative (bw derivative)  ////////////////////////
  
  if(paral_dir[0])
    {
      //split fw T border: VN 1 goes to fw out border, VN 0 goes to VN 1
      NISSA_PARALLEL_LOOP(base_src,0,bgqlx_t_vbord_vol/2)
	{
	  //the source starts at the middle of result border buffer
	  int isrc=base_src+bgqlx_t_vbord_vol/2;
	  //local shuffling: inside the result border buffer, data is stored as in the border, so base_src points
	  //already to the correct destination, but we must shift for volume and 4 terms of bw der
	  int idst_loc=8*base_src+4;
	  //non-local shuffling: must enter fw buffer for direction 0
	  int idst_buf=base_src;
	  halfspincolor temp;
	  BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(temp,((halfspincolor*)nissa_send_buf)[idst_buf],
					    bgq_hopping_matrix_output_T_buffer[isrc]);
	  HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_binded[idst_loc],temp,1); //copy VN=0 to VN=1
	}
    }
  else
    //we have only to transpose between VN
    NISSA_PARALLEL_LOOP(base_src,0,bgqlx_t_vbord_vol/2)
      {
	//the source starts at the middle of result border buffer
	int isrc=base_src+bgqlx_t_vbord_vol/2;
	int idst=8*base_src+4; //look at previous VN=0 for explanation
	BI_HALFSPINCOLOR_TRANSPOSE(bgq_hopping_matrix_output_binded[idst],bgq_hopping_matrix_output_T_buffer[isrc]);
      }
  
  //all the buffers are filled, communications can be start
  thread_barrier(HOPPING_MATRIX_APPLICATION_BARRIER);
}

//perform communications between VN and start all the communications between nodes
void start_bgq_Wilson_hopping_matrix_communications()
{
  //shuffe data between virtual nodes and fill T out buffer
  bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill();
  
  /*debug:
  GET_THREAD_ID();
  //check that all sending buffer is filled
  NISSA_PARALLEL_LOOP(i,0,bgqlx_vbord_vol*2*3*2)
  if(((complex*)nissa_send_buf)[i]) crash("");*/
  
  //start buffered communications of scattered data to other nodes
  buffered_comm_start(buffered_lx_spincolor_comm);
}

//finish the communications and put in place the communicated data
void finish_bgq_Wilson_hopping_matrix_communications()
{
  GET_THREAD_ID();
  
  //wait the end of communications
  buffered_comm_wait(buffered_lx_spincolor_comm);
  
  /*debug:
  //check that all sending buffer is filled
  NISSA_PARALLEL_LOOP(i,0,bgqlx_vbord_vol*2*3*2)
  if(((complex*)nissa_recv_buf)[i]) crash("");*/

  //T bw border (bw derivative): data goes to VN 0
  NISSA_PARLLEL_LOOP(isrc,0,bgqlx_t_vbord_vol/2)
    {
      int idst=
      HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(bgq_hopping_matrix_output_binded[idst],((halfspincolor*)nissa_recv_buf)[isrc],0);
    }
  
  //other 3 bw borders

  //T fw border (fw derivative)
  
  //other 3 fw borders
}
