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

#define HOP_HEADER(A)				\
  REORDER_BARRIER();				\
  CACHE_PREFETCH(out+A);			\
  BI_SU3_PREFETCH_NEXT(links[A])

#define LOAD_AND_SUMM_NEXT_TERM(OUT,TEMP,A)	\
  REORDER_BARRIER();				\
  REG_LOAD_BI_COLOR_ADVANCING(TEMP,A);		\
  REG_BI_COLOR_SUMM(OUT,OUT,TEMP);

#define LOAD_AND_SUBT_NEXT_TERM(OUT,TEMP,A)	\
  REORDER_BARRIER();				\
  REG_LOAD_BI_COLOR_ADVANCING(TEMP,A);		\
  REG_BI_COLOR_SUBT(OUT,OUT,TEMP);

#ifdef BGQ
 #define SITE_COPY(base_out,base_in)		\
  {						\
    void *in=base_in,*out=base_out;		\
    DECLARE_REG_BI_COLOR(reg_temp);		\
    BI_COLOR_PREFETCH_NEXT(in);			\
    REG_LOAD_BI_COLOR(reg_temp,in);		\
    STORE_REG_BI_COLOR(out,reg_temp);		\
  }
#else
 #define SITE_COPY(out,in) BI_COLOR_COPY(out,in)
#endif

//summ the eight contributions and divide by two
THREADABLE_FUNCTION_1ARG(hopping_matrix_eo_or_eo_expand_to_D, bi_color*,out)
{
  GET_THREAD_ID();
  
  //wait that all the terms are put in place
  THREAD_BARRIER();
  
  //result of split application
  bi_color *bgq_hopping_matrix_output_data=(bi_color*)nissa_send_buf+bord_volh/2;
  
  //define workload and point to the begin of each chunk
  NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,loc_volh/2,thread_id,NACTIVE_THREADS);
  void *temp_ptr=(bi_complex*)(bgq_hopping_matrix_output_data+start*8)-1;
  void *out_ptr=(bi_complex*)(out+start)-1;
  
  //regs
  DECLARE_REG_BI_COLOR(reg_out);
  DECLARE_REG_BI_COLOR(reg_temp);
  DECLARE_REG_BI_COMPLEX(reg_one_half);
  REG_SPLAT_BI_COMPLEX(reg_one_half,0.5);
  
  for(int i=start;i<end;i++)
    {
      //copy first term
      REG_LOAD_BI_COLOR_ADVANCING(reg_out,temp_ptr);
      
      //other 7 terms
      LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	
      //put final 0.5 (not with the minus!)
      REG_BI_COLOR_PROD_4DOUBLE(reg_out,reg_out,reg_one_half);
      STORE_REG_BI_COLOR_ADVANCING(out_ptr,reg_out);
    }

  //final sync
  set_borders_invalid(out);
}}

//summ the eight contributions, divide by two and subtract from the diagonal squared mass term
THREADABLE_FUNCTION_3ARG(hopping_matrix_eo_or_eo_expand_to_D_subtract_from_mass2_times_in, bi_color*,out, double,mass2, bi_color*,in)
{
  GET_THREAD_ID();
  
  //wait that all the terms are put in place
  THREAD_BARRIER();
  
  //result of split application
  bi_color *bgq_hopping_matrix_output_data=(bi_color*)nissa_send_buf+bord_volh/2;
  
  //define workload and point to the begin of each chunk
  NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,loc_volh/2,thread_id,NACTIVE_THREADS);
  void *temp_ptr=(bi_complex*)(bgq_hopping_matrix_output_data+start*8)-1;
  void *out_ptr=(bi_complex*)(out+start)-1;
  void *in_ptr=(bi_complex*)(in+start)-1;
  
  //regs
  DECLARE_REG_BI_COLOR(reg_in);
  DECLARE_REG_BI_COLOR(reg_out);
  DECLARE_REG_BI_COLOR(reg_temp);
  
  //-0.5
  DECLARE_REG_BI_COMPLEX(reg_mone_half);
  REG_SPLAT_BI_COMPLEX(reg_mone_half,-0.5);
  
  //reg_mass2
  DECLARE_REG_BI_COMPLEX(reg_mass2);
  REG_SPLAT_BI_COMPLEX(reg_mass2,mass2);
  
  for(int i=start;i<end;i++)
    {
      //copy first term
      REG_LOAD_BI_COLOR_ADVANCING(reg_out,temp_ptr);
      
      //other 7 terms
      LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
      LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	
      //put final -0.5
      REG_BI_COLOR_PROD_4DOUBLE(reg_out,reg_out,reg_mone_half);
      
      //load diagonal term, and summ it
      REG_LOAD_BI_COLOR_ADVANCING(reg_in,in_ptr);
      REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(reg_out,reg_out,reg_in,reg_mass2);
      
      //store
      STORE_REG_BI_COLOR_ADVANCING(out_ptr,reg_out);
    }

  //final sync
  set_borders_invalid(out);
}}

THREADABLE_FUNCTION_5ARG(apply_staggered_hopping_matrix_oe_or_eo_bgq_nocomm_nobarrier, bi_oct_su3**,conf, int,istart, int,iend, bi_color*,in, int,oe_or_eo)
{
  GET_THREAD_ID();
  
  bi_color *out=(bi_color*)nissa_send_buf;
  
  NISSA_PARALLEL_LOOP(ibgqlx,istart,iend)
    {
      //take short access to link and output indexing
      int *iout=viroe_or_vireo_hopping_matrix_output_pos[oe_or_eo].inter_fr_in_pos+ibgqlx*8;
      bi_su3 *links=(bi_su3*)(conf[oe_or_eo]+ibgqlx);
      
      //declare and load
      DECLARE_REG_BI_COLOR(reg_in);
      REG_LOAD_BI_COLOR(reg_in,in[ibgqlx]);
      
      HOP_HEADER(0); //T backward scatter (forward derivative)
      REG_BI_SU3_PROD_BI_COLOR_LOAD_STORE(out[iout[0]],links[0],reg_in);
      HOP_HEADER(1); //X backward scatter (forward derivative)
      REG_BI_SU3_PROD_BI_COLOR_LOAD_STORE(out[iout[1]],links[1],reg_in);
      HOP_HEADER(2); //Y backward scatter (forward derivative)
      REG_BI_SU3_PROD_BI_COLOR_LOAD_STORE(out[iout[2]],links[2],reg_in);
      HOP_HEADER(3); //Z backward scatter (forward derivative)
      REG_BI_SU3_PROD_BI_COLOR_LOAD_STORE(out[iout[3]],links[3],reg_in);
      
      HOP_HEADER(4); //T forward scatter (backward derivative)
      REG_BI_SU3_DAG_PROD_BI_COLOR_LOAD_STORE(out[iout[4]],links[4],reg_in);
      HOP_HEADER(5); //X forward scatter (backward derivative)
      REG_BI_SU3_DAG_PROD_BI_COLOR_LOAD_STORE(out[iout[5]],links[5],reg_in);
      HOP_HEADER(6); //Y forward scatter (backward derivative)
      REG_BI_SU3_DAG_PROD_BI_COLOR_LOAD_STORE(out[iout[6]],links[6],reg_in);
      HOP_HEADER(7); //Z forward scatter (backward derivative)
      REG_BI_SU3_DAG_PROD_BI_COLOR_LOAD_STORE(out[iout[7]],links[7],reg_in);
    }
}}

//if virtual parallelized dir is really parallelized, fill send buffers
THREADABLE_FUNCTION_0ARG(bgq_staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill)
{
  GET_THREAD_ID();
  
  //vdir buffers might have been filled in another order
  THREAD_BARRIER();
  
  //form two thread team
  FORM_TWO_THREAD_TEAMS();
  
  //short access
  int v=nissa_vnode_paral_dir;
  const int fact=2;
  bi_color *bgq_hopping_matrix_output_data=(bi_color*)nissa_send_buf+bord_volh/fact; //half vol bisp = vol sp 
  bi_color *bgq_hopping_matrix_output_vdir_buffer=bgq_hopping_matrix_output_data+8*loc_volh/fact;
  
  ///////////////////////// bw scattered v derivative (fw derivative)  ////////////////////////
  
  if(is_in_first_team)
    //split bw v border: VN 0 goes to bw out border (first half)
    NISSA_CHUNK_LOOP(base_isrc,0,vbord_vol/4/fact,thread_in_team_id,nthreads_in_team)
      {
	//the source starts at the middle of result border buffer
	int isrc=2*base_isrc+0*vbord_vol/2/fact; //we match 2 sites
	//non-local shuffling: must enter bw buffer for direction v
	int idst_buf=(0*bord_volh/2+bord_offset[v]/2)/fact+base_isrc;
	//load the first
	DECLARE_REG_BI_COLOR(in0);
	REG_LOAD_BI_COLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	//load the second
	DECLARE_REG_BI_COLOR(in1);
	REG_LOAD_BI_COLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	//merge the two and save
	DECLARE_REG_BI_COLOR(to_buf);
	REG_BI_COLOR_V0_MERGE(to_buf,in0,in1);
	STORE_REG_BI_COLOR(((bi_color*)nissa_send_buf)[idst_buf],to_buf);
      }
  
  ///////////////////////// fw scattered v derivative (bw derivative)  ////////////////////////
  
  if(is_in_second_team)
    //split fw v border: VN 1 goes to fw out border (second half)
    NISSA_CHUNK_LOOP(base_isrc,0,vbord_vol/4/fact,thread_in_team_id,nthreads_in_team)
      {
	//the source starts at the middle of result border buffer
	int isrc=2*base_isrc+1*vbord_vol/2/fact;
	//non-local shuffling: must enter fw buffer (starting at bord_volh/2 because its bi) for direction v
	int idst_buf=(1*bord_volh/2+bord_offset[v]/2)/fact+base_isrc;
	//load the first
	DECLARE_REG_BI_COLOR(in0);
	REG_LOAD_BI_COLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	//load the second
	DECLARE_REG_BI_COLOR(in1);
	REG_LOAD_BI_COLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	//merge the two and save
	DECLARE_REG_BI_COLOR(to_buf);
	REG_BI_COLOR_V1_MERGE(to_buf,in0,in1);
	STORE_REG_BI_COLOR(((bi_color*)nissa_send_buf)[idst_buf],to_buf);
      }
}}

//pick data from vbuffer and put it in correct position
THREADABLE_FUNCTION_0ARG(bgq_staggered_hopping_matrix_oe_or_eo_vdir_VN_local_transpose)
{
  GET_THREAD_ID();
  
  //short access
  const int fact=2;
  int v=nissa_vnode_paral_dir;
  bi_color *bgq_hopping_matrix_output_vdir_buffer=(bi_color*)nissa_send_buf+bord_volh/fact+8*loc_volh/fact;
  
  FORM_TWO_THREAD_TEAMS();

  //backward scatter (forward derivative)
  bi_color *base_out_bw=(bi_color*)nissa_send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
  bi_color *base_in_bw=bgq_hopping_matrix_output_vdir_buffer+0*vdir_bord_vol/fact;
  if(is_in_first_team)
    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
      {
	//load
	DECLARE_REG_BI_COLOR(reg_in);
	REG_LOAD_BI_COLOR(reg_in,base_in_bw[isrc]);
	//transpose and store
	DECLARE_REG_BI_COLOR(reg_out);
	REG_BI_COLOR_TRANSPOSE(reg_out,reg_in);
	STORE_REG_BI_COLOR(base_out_bw[isrc*8+0+v],reg_out);
      }
  
  //forward scatter (backward derivative)
  bi_color *base_out_fw=(bi_color*)nissa_send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;  
  bi_color *base_in_fw=bgq_hopping_matrix_output_vdir_buffer+1*vdir_bord_vol/fact;
  if(is_in_second_team)
    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
      {
	//load
	DECLARE_REG_BI_COLOR(reg_in);
	REG_LOAD_BI_COLOR(reg_in,base_in_fw[isrc]);
	//transpose and store
	DECLARE_REG_BI_COLOR(reg_out);
	REG_BI_COLOR_TRANSPOSE(reg_out,reg_in);
	STORE_REG_BI_COLOR(base_out_fw[isrc*8+4+v],reg_out);
      }
}}
  
//perform communications between VN and start all the communications between nodes
THREADABLE_FUNCTION_0ARG(start_staggered_hopping_matrix_oe_or_eo_bgq_communications)
{
  //shuffle data between virtual nodes and fill vdir out buffer
  if(paral_dir[nissa_vnode_paral_dir]) bgq_staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill();
  
  //after the barrier, all buffers are filled and communications can start
  THREAD_BARRIER();
  
  //start communications of scattered data to other nodes
  comm_start(eo_color_comm);
  
  //if v dir is not parallelized we have only to transpose between VN, and no real communication happens
  if(!paral_dir[nissa_vnode_paral_dir]) bgq_staggered_hopping_matrix_oe_or_eo_vdir_VN_local_transpose();
}}

//finish the communications and put in place the communicated data
THREADABLE_FUNCTION_1ARG(finish_staggered_hopping_matrix_oe_or_eo_bgq_communications, int,oe_or_eo)
{
  GET_THREAD_ID();
  FORM_TWO_THREAD_TEAMS();
  
  //short access
  int v=nissa_vnode_paral_dir;
  const int fact=2;
  int *def_pos=viroe_or_vireo_hopping_matrix_output_pos[oe_or_eo].inter_fr_recv_pos;
  
  //wait communications end
  comm_wait(eo_color_comm);
  
  //vdir bw border (bw derivative): data goes to VN 0
  if(is_in_first_team)
    if(paral_dir[v])
      {
	//inside incoming borders vdir is ordered naturally, while in the output data it comes first
	bi_color *base_out=(bi_color*)nissa_send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;
	bi_color *base_vdir_in=(bi_color*)nissa_send_buf+(bord_volh+8*loc_volh+1*vbord_vol/2)/fact;
	bi_color *base_bord_in=(bi_color*)nissa_recv_buf+bord_offset[v]/2/fact;
	
	NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/2/fact,thread_in_team_id,nthreads_in_team)
	  {
	    //VN=0 must be filled with border
	    DECLARE_REG_BI_COLOR(in0);
	    REG_LOAD_BI_COLOR(in0,base_bord_in[isrc]);
	    //VN=1 with buf0
	    DECLARE_REG_BI_COLOR(in1);
	    REG_LOAD_BI_COLOR(in1,base_vdir_in[2*isrc]);
	    //merge and save
	    DECLARE_REG_BI_COLOR(to_dest);
	    REG_BI_COLOR_V0_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_COLOR(base_out[2*isrc*8+4+v],to_dest);
	    
	    //VN=1 with buf1
	    REG_LOAD_BI_COLOR(in1,base_vdir_in[2*isrc+1]);
	    //merge and save
	    REG_BI_COLOR_V10_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_COLOR(base_out[(2*isrc+1)*8+4+v],to_dest);
	  }
      }
  
  //other 3 bw borders
  if(is_in_first_team)
    for(int imu=0;imu<3;imu++)
      {
	int mu=perp_dir[v][imu];
	bi_color *base_out=(bi_color*)nissa_send_buf+bord_volh/fact;
	bi_color *base_in=(bi_color*)nissa_recv_buf;
	NISSA_CHUNK_LOOP(isrc,bord_offset[mu]/2/fact,(bord_offset[mu]+bord_dir_vol[mu])/2/fact,
			 thread_in_team_id,nthreads_in_team)
	  SITE_COPY(base_out[def_pos[isrc]],base_in[isrc]);
      }
  
  //v fw border (fw derivative): data goes to VN 1
  if(is_in_second_team)
    if(paral_dir[v])
      {
	//inside incoming borders vdir is ordered naturally, while in the output data it comes first
	bi_color *base_out=(bi_color*)nissa_send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
	bi_color *base_vdir_in=(bi_color*)nissa_send_buf+(bord_volh+8*loc_volh+0*vbord_vol/2)/fact;
	bi_color *base_bord_in=(bi_color*)nissa_recv_buf+(bord_vol/4+bord_offset[v]/2)/fact;
	
	NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/2/fact,thread_in_team_id,nthreads_in_team)
	  {
	    //VN=0 with buf1
	    DECLARE_REG_BI_COLOR(in0);
	    REG_LOAD_BI_COLOR(in0,base_vdir_in[2*isrc]);
	    //VN=1 must be filled with border 0
	    DECLARE_REG_BI_COLOR(in1);
	    REG_LOAD_BI_COLOR(in1,base_bord_in[isrc]);
	    //merge and save
	    DECLARE_REG_BI_COLOR(to_dest);
	    REG_BI_COLOR_V10_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_COLOR(base_out[2*isrc*8+0+v],to_dest);
	    
	    //VN=0 with buf1
	    REG_LOAD_BI_COLOR(in0,base_vdir_in[2*isrc+1]);
	    //merge and save
	    REG_BI_COLOR_V1_MERGE(to_dest,in0,in1);
	    STORE_REG_BI_COLOR(base_out[(2*isrc+1)*8+0+v],to_dest);
	  }
      }
  
  //other 3 fw borders
  if(is_in_second_team)
    for(int imu=0;imu<3;imu++)
      {
	int mu=perp_dir[v][imu];
	bi_color *base_out=(bi_color*)nissa_send_buf+bord_volh/fact;
	bi_color *base_in=(bi_color*)nissa_recv_buf;
	NISSA_CHUNK_LOOP(isrc,(bord_volh+bord_offset[mu])/2/fact,(bord_volh+bord_offset[mu]+bord_dir_vol[mu])/2/fact,
			 thread_in_team_id,nthreads_in_team)
	  SITE_COPY(base_out[def_pos[isrc]],base_in[isrc]);
      }
}}
