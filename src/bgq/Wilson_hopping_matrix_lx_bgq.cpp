#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "communicate/borders.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "bgq_macros.hpp"

/*
  In bgq version we merge two sites along t directions, that is, (t,x,y,z) and (t+T/2,x,y,z),
  so that only site with time coordinate between 0 and T/2-1 must be considered.
  Since we want to load everything sequentially, we need to duplicate the gauge configuration.

  We apply hopping matrix scanning on sink index, then store the 8 contributions to each source separately.
  They can be summed outside according to external usage.
  First 4 entries contains forward derivative scattered backward, then others.
  
  We do not sync and do not perform any communication, so this must be done outside.
*/

#define PROJ_HEADER(A)				\
  REORDER_BARRIER();				\
  CACHE_PREFETCH(out+A);			\
  BI_SU3_PREFETCH_NEXT(links[A])

#ifdef BGQ
 #define SITE_COPY(base_out,base_in)			\
  {							\
    void *in=base_in,*out=base_out;			\
    DECLARE_REG_BI_HALFSPINCOLOR(reg_temp);		\
    BI_HALFSPINCOLOR_PREFETCH_NEXT(in);			\
    REG_LOAD_BI_HALFSPINCOLOR(reg_temp,in);		\
    STORE_REG_BI_HALFSPINCOLOR(out,reg_temp);		\
  }
#else
 #define SITE_COPY(out,in) BI_HALFSPINCOLOR_COPY(out,in)
#endif

THREADABLE_FUNCTION_4ARG(apply_Wilson_hopping_matrix_lx_bgq_nocomm_nobarrier, bi_oct_su3*,conf, int,istart, int,iend, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  bi_halfspincolor *out=(bi_halfspincolor*)nissa_send_buf;
  
  NISSA_PARALLEL_LOOP(ibgqlx,istart,iend)
    {
      //take short access to link and output indexing
      int *iout=virlx_hopping_matrix_output_pos.inter_fr_in_pos+ibgqlx*8;
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
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[0]],links[0],reg_proj);
      
      //X backward scatter (forward derivative)
      PROJ_HEADER(1);
      REG_BI_COLOR_ISUMM(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_ISUMM(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[1]],links[1],reg_proj);
      
      //Y backward scatter (forward derivative)
      PROJ_HEADER(2);
      REG_BI_COLOR_SUMM(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_SUBT(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[2]],links[2],reg_proj);
      
      //Z backward scatter (forward derivative)
      PROJ_HEADER(3);
      REG_BI_COLOR_ISUMM(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_ISUBT(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[3]],links[3],reg_proj);
      
      //T forward scatter (backward derivative)
      PROJ_HEADER(4);
      REG_BI_COLOR_SUBT(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_SUBT(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[4]],links[4],reg_proj);
      
      //X forward scatter (backward derivative)
      PROJ_HEADER(5);
      REG_BI_COLOR_ISUBT(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_ISUBT(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[5]],links[5],reg_proj);

      //Y forward scatter (backward derivative)
      PROJ_HEADER(6);
      REG_BI_COLOR_SUBT(reg_proj_s0,reg_in_s0,reg_in_s3);
      REG_BI_COLOR_SUMM(reg_proj_s1,reg_in_s1,reg_in_s2);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[6]],links[6],reg_proj);
      
      //Z forward scatter (backward derivative)
      PROJ_HEADER(7);
      REG_BI_COLOR_ISUBT(reg_proj_s0,reg_in_s0,reg_in_s2);
      REG_BI_COLOR_ISUMM(reg_proj_s1,reg_in_s1,reg_in_s3);
      REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out[iout[7]],links[7],reg_proj);
    }
}}
  
//swap border between VN and, if virtual parallelized dir is really parallelized, fill send buffers
THREADABLE_FUNCTION_0ARG(bgq_Wilson_hopping_matrix_lx_vdir_VN_comm_and_buff_fill)
{
  GET_THREAD_ID();
  
  //vdir buffers might have been filled in another order
  THREAD_BARRIER();
  
  //form two thread team
  FORM_TWO_THREAD_TEAMS();
  
  //short access
  int v=nissa_vnode_paral_dir;
  bi_halfspincolor *bgq_hopping_matrix_output_data=(bi_halfspincolor*)nissa_send_buf+bord_volh; //half vol bisp = vol sp 
  bi_halfspincolor *bgq_hopping_matrix_output_vdir_buffer=bgq_hopping_matrix_output_data+8*loc_volh;
  
  ///////////////////////// bw scattered v derivative (fw derivative)  ////////////////////////
  
  if(is_in_first_team)
    //split bw v border: VN 0 goes to bw out border (first half)
    NISSA_CHUNK_LOOP(base_isrc,0,vbord_vol/4,thread_in_team_id,nthreads_in_team)
      {
	//the source starts at the middle of result border buffer
	int isrc=2*base_isrc+0*vbord_vol/2; //we match 2 sites
	//non-local shuffling: must enter bw buffer for direction v
	int idst_buf=0*bord_volh/2+bord_offset[v]/2+base_isrc;
	//load the first
	DECLARE_REG_BI_HALFSPINCOLOR(in0);
	REG_LOAD_BI_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	//load the second
	DECLARE_REG_BI_HALFSPINCOLOR(in1);
	REG_LOAD_BI_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	//merge the two and save
	DECLARE_REG_BI_HALFSPINCOLOR(to_buf);
	REG_BI_HALFSPINCOLOR_V0_MERGE(to_buf,in0,in1);
	STORE_REG_BI_HALFSPINCOLOR(((bi_halfspincolor*)nissa_send_buf)[idst_buf],to_buf);
      }

  ///////////////////////// fw scattered v derivative (bw derivative)  ////////////////////////
  
  if(is_in_second_team)
    //split fw v border: VN 1 goes to fw out border (second half)
    NISSA_CHUNK_LOOP(base_isrc,0,vbord_vol/4,thread_in_team_id,nthreads_in_team)
      {
	//the source starts at the middle of result border buffer
	int isrc=2*base_isrc+1*vbord_vol/2;
	//non-local shuffling: must enter fw buffer (starting at bord_volh/2 because its bi) for direction 0
	int idst_buf=1*bord_volh/2+bord_offset[v]/2+base_isrc;
	//load the first
	DECLARE_REG_BI_HALFSPINCOLOR(in0);
	REG_LOAD_BI_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	//load the second
	DECLARE_REG_BI_HALFSPINCOLOR(in1);
	REG_LOAD_BI_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	//merge the two and save
	DECLARE_REG_BI_HALFSPINCOLOR(to_buf);
	REG_BI_HALFSPINCOLOR_V1_MERGE(to_buf,in0,in1);
	STORE_REG_BI_HALFSPINCOLOR(((bi_halfspincolor*)nissa_send_buf)[idst_buf],to_buf);
      }
}}

//pick data from vbuffer and put it in correct position
THREADABLE_FUNCTION_0ARG(bgq_Wilson_hopping_matrix_lx_vdir_VN_transpose)
{
  GET_THREAD_ID();
  
  //short access
  const int fact=1;
  int v=nissa_vnode_paral_dir;
  bi_halfspincolor *bgq_hopping_matrix_output_vdir_buffer=(bi_halfspincolor*)nissa_send_buf+bord_volh/fact+
    8*loc_volh/fact;
  
  FORM_TWO_THREAD_TEAMS();

  //backward scatter (forward derivative)
  bi_halfspincolor *base_out_bw=(bi_halfspincolor*)nissa_send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
  bi_halfspincolor *base_in_bw=bgq_hopping_matrix_output_vdir_buffer+0*vdir_bord_vol/fact;
  if(is_in_first_team)
    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
      {
        //load
        DECLARE_REG_BI_HALFSPINCOLOR(reg_in);
        REG_LOAD_BI_HALFSPINCOLOR(reg_in,base_in_bw[isrc]);
        //transpose and store
        DECLARE_REG_BI_HALFSPINCOLOR(reg_out);
        REG_BI_HALFSPINCOLOR_TRANSPOSE(reg_out,reg_in);
        STORE_REG_BI_HALFSPINCOLOR(base_out_bw[isrc*8+0+v],reg_out);
      }
  
  //forward scatter (backward derivative)
  bi_halfspincolor *base_out_fw=(bi_halfspincolor*)nissa_send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;  
  bi_halfspincolor *base_in_fw=bgq_hopping_matrix_output_vdir_buffer+1*vdir_bord_vol/fact;
  if(is_in_second_team)
    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
      {
        //load
        DECLARE_REG_BI_HALFSPINCOLOR(reg_in);
        REG_LOAD_BI_HALFSPINCOLOR(reg_in,base_in_fw[isrc]);
        //transpose and store
        DECLARE_REG_BI_HALFSPINCOLOR(reg_out);
        REG_BI_HALFSPINCOLOR_TRANSPOSE(reg_out,reg_in);
        STORE_REG_BI_HALFSPINCOLOR(base_out_fw[isrc*8+4+v],reg_out);
      }
}}

//perform communications between VN and start all the communications between nodes
THREADABLE_FUNCTION_0ARG(start_Wilson_hopping_matrix_lx_bgq_communications)
{
  //shuffle data between virtual nodes and fill vdir out buffer
  if(paral_dir[nissa_vnode_paral_dir]) bgq_Wilson_hopping_matrix_lx_vdir_VN_comm_and_buff_fill();
  
  //after the barrier, all buffers are filled and communications can start
  THREAD_BARRIER();
  
  //start communications of scattered data to other nodes
  comm_start(lx_halfspincolor_comm);
  
  //if v dir is not parallelized we have only to transpose between VN, and no real communication happens
  if(!paral_dir[nissa_vnode_paral_dir]) bgq_Wilson_hopping_matrix_lx_vdir_VN_transpose();  
}}

//finish the communications and put in place the communicated data
THREADABLE_FUNCTION_0ARG(finish_Wilson_hopping_matrix_lx_bgq_communications)
{
  GET_THREAD_ID();
  FORM_TWO_THREAD_TEAMS();
  
  //short access
  int v=nissa_vnode_paral_dir;
  int *def_pos=virlx_hopping_matrix_output_pos.inter_fr_recv_pos;
  
  //wait communications end
  comm_wait(lx_halfspincolor_comm);
  
  //vdir bw border (bw derivative): data goes to VN 0
  if(paral_dir[v])
    if(is_in_first_team)
      {
	//inside incoming borders vdir is ordered naturally, while in the output data it comes first
	bi_halfspincolor *base_out=(bi_halfspincolor*)nissa_send_buf+bord_volh+0*8*bord_dir_vol[v];
	bi_halfspincolor *base_vdir_in=(bi_halfspincolor*)nissa_send_buf+bord_volh+8*loc_volh+1*vbord_vol/2;
	bi_halfspincolor *base_bord_in=(bi_halfspincolor*)nissa_recv_buf+bord_offset[v]/2;
	
	NISSA_CHUNK_LOOP(isrc,0,bord_dir_vol[v]/2,thread_in_team_id,nthreads_in_team)
	  {
	    //VN=0 must be filled with border
	    DECLARE_REG_BI_HALFSPINCOLOR(in0);
	    REG_LOAD_BI_HALFSPINCOLOR(in0,base_bord_in[isrc]);
	    //VN=1 with buf0
	    DECLARE_REG_BI_HALFSPINCOLOR(in1);
	    REG_LOAD_BI_HALFSPINCOLOR(in1,base_vdir_in[2*isrc]);
	    //merge and save
	    DECLARE_REG_BI_HALFSPINCOLOR(to_dest);
	    REG_BI_HALFSPINCOLOR_V0_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_HALFSPINCOLOR(base_out[2*isrc*8+4+v],to_dest);
	    
	    //VN=1 with buf1
	    REG_LOAD_BI_HALFSPINCOLOR(in1,base_vdir_in[2*isrc+1]);
	    //merge and save
	    REG_BI_HALFSPINCOLOR_V10_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_HALFSPINCOLOR(base_out[(2*isrc+1)*8+4+v],to_dest);
	  }
      }
  
  //other 3 bw borders
  if(is_in_first_team)
    for(int imu=0;imu<3;imu++)
      {
	int mu=perp_dir[v][imu];
	bi_halfspincolor *base_out=(bi_halfspincolor*)nissa_send_buf+bord_volh;
	bi_halfspincolor *base_in=(bi_halfspincolor*)nissa_recv_buf;
	NISSA_CHUNK_LOOP(isrc,bord_offset[mu]/2,bord_offset[mu]/2+bord_dir_vol[mu]/2,
			 thread_in_team_id,nthreads_in_team)
	  SITE_COPY(base_out[def_pos[isrc]],base_in[isrc]);
      }
  
  //v fw border (fw derivative): data goes to VN 1
  if(is_in_second_team)
    if(paral_dir[v])
      {
	//inside incoming borders vdir is ordered naturally, while in the output data it comes first
	bi_halfspincolor *base_out=(bi_halfspincolor*)nissa_send_buf+bord_volh+1*8*bord_dir_vol[v];
	bi_halfspincolor *base_vdir_in=(bi_halfspincolor*)nissa_send_buf+bord_volh+8*loc_volh+0*vbord_vol/2;
	bi_halfspincolor *base_bord_in=(bi_halfspincolor*)nissa_recv_buf+bord_vol/4+bord_offset[v]/2;
	
	NISSA_CHUNK_LOOP(isrc,0,bord_dir_vol[v]/2,thread_in_team_id,nthreads_in_team)
	  {
	    //VN=0 with buf1
	    DECLARE_REG_BI_HALFSPINCOLOR(in0);
	    REG_LOAD_BI_HALFSPINCOLOR(in0,base_vdir_in[2*isrc]);
	    //VN=1 must be filled with border 0
	    DECLARE_REG_BI_HALFSPINCOLOR(in1);
	    REG_LOAD_BI_HALFSPINCOLOR(in1,base_bord_in[isrc]);
	    //merge and save
	    DECLARE_REG_BI_HALFSPINCOLOR(to_dest);
	    REG_BI_HALFSPINCOLOR_V10_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_HALFSPINCOLOR(base_out[2*isrc*8+0+v],to_dest);
	    
	    //VN=0 with buf1
	    REG_LOAD_BI_HALFSPINCOLOR(in0,base_vdir_in[2*isrc+1]);
	    //merge and save
	    REG_BI_HALFSPINCOLOR_V1_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_HALFSPINCOLOR(base_out[(2*isrc+1)*8+0+v],to_dest);
	}
      }
  
  //other 3 fw borders
  if(is_in_second_team)
    for(int imu=0;imu<3;imu++)
      {
	int mu=perp_dir[v][imu];
	bi_halfspincolor *base_out=(bi_halfspincolor*)nissa_send_buf+bord_volh;
	bi_halfspincolor *base_in=(bi_halfspincolor*)nissa_recv_buf;
	NISSA_CHUNK_LOOP(isrc,bord_vol/4+bord_offset[mu]/2,bord_vol/4+bord_offset[mu]/2+bord_dir_vol[mu]/2,
			 thread_in_team_id,nthreads_in_team)
	  SITE_COPY(base_out[def_pos[isrc]],base_in[isrc]);
      }
}}
