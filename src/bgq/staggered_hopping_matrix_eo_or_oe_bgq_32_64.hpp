#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_vir.hpp"
#include "new_types/complex.hpp"
#include "Wilson_hopping_matrix_lx_bgq.hpp"
#include "threads/threads.hpp"

#include "bgq_macros.hpp"

#if (((defined VERSION_64) && (defined VERSION_32)) || ((!defined VERSION_64) && (!defined VERSION_32)))
 #error Impossible!
#endif

#if ((defined VERSION_64) && (!defined VERSION_32))
 #define PREC_TYPE double
 #define VIR_32_64_COMPLEX vir_complex
 #define VIR_32_64_COLOR vir_color
 #define VIR_32_64_SU3 vir_su3
 #define VIR_32_64_OCT_SU3 vir_oct_su3
 #define EO_32_64_COLOR_COMM eo_color_comm
 #define BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL \
         bgq_double_staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill
 #define BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE \
         bgq_double_staggered_hopping_matrix_oe_or_eo_vdir_VN_local_transpose
 #define START_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         start_double_staggered_hopping_matrix_oe_or_eo_bgq_communications
 #define FINISH_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         finish_double_staggered_hopping_matrix_oe_or_eo_bgq_communications
 #define HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D \
         hopping_matrix_oe_or_eo_expand_to_double_staggered_D_bgq
 #define HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D_SUBTRACT_FROM_MASS2_TIMES_IN \
         hopping_matrix_oe_or_eo_expand_to_double_staggered_D_subtract_from_mass2_times_in_bgq
#define APPLY_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM	\
         apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm
 #define VIR_32_64_SU3_PREFETCH_NEXT VIR_SU3_PREFETCH_NEXT
 #define REG_LOAD_VIR_32_64_COLOR_ADVANCING REG_LOAD_VIR_COLOR_ADVANCING
 #define REG_LOAD_VIR_32_64_COLOR REG_LOAD_VIR_COLOR
 #define STORE_REG_VIR_32_64_COLOR STORE_REG_VIR_COLOR
 #define STORE_REG_VIR_32_64_COLOR_ADVANCING STORE_REG_VIR_COLOR_ADVANCING
 #define VIR_32_64_COLOR_PREFETCH_NEXT VIR_COLOR_PREFETCH_NEXT
 #define REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE REG_VIR_SU3_PROD_VIR_COLOR_LOAD_STORE
 #define REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE REG_VIR_SU3_DAG_PROD_VIR_COLOR_LOAD_STORE
#endif

#if ((defined VERSION_32) && (!defined VERSION_64))
 #define PREC_TYPE float
 #define VIR_32_64_COMPLEX vir_single_complex
 #define VIR_32_64_COLOR vir_single_color
 #define VIR_32_64_SU3 vir_single_su3
 #define VIR_32_64_OCT_SU3 vir_single_oct_su3
 #define EO_32_64_COLOR_COMM eo_single_color_comm
 #define BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL \
         bgq_single_staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill
 #define BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE \
         bgq_single_staggered_hopping_matrix_oe_or_eo_vdir_VN_local_transpose
 #define START_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         start_single_staggered_hopping_matrix_oe_or_eo_bgq_communications
 #define FINISH_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         finish_single_staggered_hopping_matrix_oe_or_eo_bgq_communications
 #define HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D \
         hopping_matrix_oe_or_eo_expand_to_single_staggered_D_bgq
 #define HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D_SUBTRACT_FROM_MASS2_TIMES_IN \
         hopping_matrix_oe_or_eo_expand_to_single_staggered_D_subtract_from_mass2_times_in_bgq
#define APPLY_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM	\
         apply_single_staggered_hopping_matrix_oe_or_eo_bgq_nocomm
 #define VIR_32_64_SU3_PREFETCH_NEXT VIR_SINGLE_SU3_PREFETCH_NEXT
 #define REG_LOAD_VIR_32_64_COLOR_ADVANCING REG_LOAD_VIR_SINGLE_COLOR_ADVANCING
 #define REG_LOAD_VIR_32_64_COLOR REG_LOAD_VIR_SINGLE_COLOR
 #define STORE_REG_VIR_32_64_COLOR STORE_REG_VIR_SINGLE_COLOR
 #define STORE_REG_VIR_32_64_COLOR_ADVANCING STORE_REG_VIR_SINGLE_COLOR_ADVANCING
 #define VIR_32_64_COLOR_PREFETCH_NEXT VIR_SINGLE_COLOR_PREFETCH_NEXT
 #define REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE REG_VIR_SINGLE_SU3_PROD_VIR_SINGLE_COLOR_LOAD_STORE
 #define REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE REG_VIR_SINGLE_SU3_DAG_PROD_VIR_SINGLE_COLOR_LOAD_STORE
#endif

#define HOP_HEADER(A)				\
  REORDER_BARRIER();				\
  CACHE_PREFETCH(out+A);			\
  VIR_32_64_SU3_PREFETCH_NEXT(links[A])

#define LOAD_AND_SUMM_NEXT_TERM(OUT,TEMP,A)	\
  REORDER_BARRIER();				\
  REG_LOAD_VIR_32_64_COLOR_ADVANCING(TEMP,A);	\
  REG_VIR_COLOR_SUMM(OUT,OUT,TEMP);

#define LOAD_AND_SUBT_NEXT_TERM(OUT,TEMP,A)	\
  REORDER_BARRIER();				\
  REG_LOAD_VIR_32_64_COLOR_ADVANCING(TEMP,A);	\
  REG_VIR_COLOR_SUBT(OUT,OUT,TEMP);

#ifdef BGQ
 #define SITE_COPY(base_out,base_in)		\
  {						\
    void *in=base_in,*out=base_out;		\
    DECLARE_REG_VIR_COLOR(reg_temp);		\
    VIR_32_64_COLOR_PREFETCH_NEXT(in);		\
    REG_LOAD_VIR_32_64_COLOR(reg_temp,in);	\
    STORE_REG_VIR_32_64_COLOR(out,reg_temp);	\
  }
#else
 #define SITE_COPY(out,in) VIR_COLOR_COPY(out,in)
#endif

namespace nissa
{
  //summ the eight contributions and divide by two (note that the non-coeff option is used below)
  THREADABLE_FUNCTION_2ARG(HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D, VIR_32_64_COLOR*,out, double,coeff)
  {
    GET_THREAD_ID();
    
    //result of split application
    VIR_32_64_COLOR *bgq_hopping_matrix_output_data=(VIR_32_64_COLOR*)send_buf+bord_volh/2;
    
    //define workload and point to the begin of each chunk
    NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,loc_volh/2,thread_id,NACTIVE_THREADS);
    void *temp_ptr=(VIR_32_64_COMPLEX*)(bgq_hopping_matrix_output_data+start*8)-1;
    void *out_ptr=(VIR_32_64_COMPLEX*)(out+start)-1;
    
    //regs
    DECLARE_REG_VIR_COLOR(reg_out);
    DECLARE_REG_VIR_COLOR(reg_temp);
    DECLARE_REG_VIR_COMPLEX(reg_one_half);
    REG_SPLAT_VIR_COMPLEX(reg_one_half,0.5*coeff);
    
    for(int i=start;i<end;i++)
      {
	//copy first term
	REG_LOAD_VIR_32_64_COLOR_ADVANCING(reg_out,temp_ptr);
	
	//other 7 terms
	LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	
	//put final 0.5
	REG_VIR_COLOR_PROD_4DOUBLE(reg_out,reg_out,reg_one_half);
	STORE_REG_VIR_32_64_COLOR_ADVANCING(out_ptr,reg_out);
      }
    
    //final sync
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  void HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D(VIR_32_64_COLOR *out)
  {HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D(out,1);}
  
  //summ the eight contributions, divide by two and subtract from the diagonal squared mass term
  THREADABLE_FUNCTION_3ARG(HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D_SUBTRACT_FROM_MASS2_TIMES_IN, VIR_32_64_COLOR*,out, PREC_TYPE,mass2, VIR_32_64_COLOR*,in)
  {
    GET_THREAD_ID();
    
    //result of split application
    VIR_32_64_COLOR *bgq_hopping_matrix_output_data=(VIR_32_64_COLOR*)send_buf+bord_volh/2;
    
    //define workload and point to the begin of each chunk
    NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,loc_volh/2,thread_id,NACTIVE_THREADS);
    void *temp_ptr=(VIR_32_64_COMPLEX*)(bgq_hopping_matrix_output_data+start*8)-1;
    void *out_ptr=(VIR_32_64_COMPLEX*)(out+start)-1;
    void *in_ptr=(VIR_32_64_COMPLEX*)(in+start)-1;
    
    //regs
    DECLARE_REG_VIR_COLOR(reg_in);
    DECLARE_REG_VIR_COLOR(reg_out);
    DECLARE_REG_VIR_COLOR(reg_temp);
    
    //-0.5
    DECLARE_REG_VIR_COMPLEX(reg_mone_half);
    REG_SPLAT_VIR_COMPLEX(reg_mone_half,-0.5);
    
    //reg_mass2
    DECLARE_REG_VIR_COMPLEX(reg_mass2);
    REG_SPLAT_VIR_COMPLEX(reg_mass2,mass2);
    
    for(int i=start;i<end;i++)
      {
	//copy first term
	REG_LOAD_VIR_32_64_COLOR_ADVANCING(reg_out,temp_ptr);
	
	//other 7 terms
	LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUMM_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	LOAD_AND_SUBT_NEXT_TERM(reg_out,reg_temp,temp_ptr);
	
	//put final -0.5
	REG_VIR_COLOR_PROD_4DOUBLE(reg_out,reg_out,reg_mone_half);
	
	//load diagonal term, and summ it
	REG_LOAD_VIR_32_64_COLOR_ADVANCING(reg_in,in_ptr);
	REG_VIR_COLOR_SUMM_THE_PROD_4DOUBLE(reg_out,reg_out,reg_in,reg_mass2);
	
	//store
	STORE_REG_VIR_32_64_COLOR_ADVANCING(out_ptr,reg_out);
      }
    
    //final sync
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_5ARG(APPLY_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM, VIR_32_64_OCT_SU3**,conf, int,istart, int,iend, VIR_32_64_COLOR*,in, int,oe_or_eo)
  {
    GET_THREAD_ID();
    
    VIR_32_64_COLOR *out=(VIR_32_64_COLOR*)send_buf;
    
    NISSA_PARALLEL_LOOP(ibgqlx,istart,iend)
      {
	//take short access to link and output indexing
	int *iout=viroe_or_vireo_hopping_matrix_output_pos[oe_or_eo].inter_fr_in_pos+ibgqlx*8;
	VIR_32_64_SU3 *links=(VIR_32_64_SU3*)(conf[oe_or_eo]+ibgqlx);
	
	//declare and load
	DECLARE_REG_VIR_COLOR(reg_in);
	REG_LOAD_VIR_32_64_COLOR(reg_in,in[ibgqlx]);
	
	HOP_HEADER(0); //T backward scatter (forward derivative)
	REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[0]],links[0],reg_in);
	HOP_HEADER(1); //X backward scatter (forward derivative)
	REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[1]],links[1],reg_in);
	HOP_HEADER(2); //Y backward scatter (forward derivative)
	REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[2]],links[2],reg_in);
	HOP_HEADER(3); //Z backward scatter (forward derivative)
	REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[3]],links[3],reg_in);
	
	HOP_HEADER(4); //T forward scatter (backward derivative)
	REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[4]],links[4],reg_in);
	HOP_HEADER(5); //X forward scatter (backward derivative)
	REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[5]],links[5],reg_in);
	HOP_HEADER(6); //Y forward scatter (backward derivative)
	REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[6]],links[6],reg_in);
	HOP_HEADER(7); //Z forward scatter (backward derivative)
	REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE(out[iout[7]],links[7],reg_in);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //if virtual parallelized dir is really parallelized, fill send buffers
  THREADABLE_FUNCTION_0ARG(BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL)
  {
    GET_THREAD_ID();
    
    //form two thread team
    FORM_TWO_THREAD_TEAMS();
    
    //short access
    int v=vnode_paral_dir;
    const int fact=2;
    VIR_32_64_COLOR *bgq_hopping_matrix_output_data=(VIR_32_64_COLOR*)send_buf+bord_volh/fact; //half vol bisp = vol sp 
    VIR_32_64_COLOR *bgq_hopping_matrix_output_vdir_buffer=bgq_hopping_matrix_output_data+8*loc_volh/fact;
    
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
	  DECLARE_REG_VIR_COLOR(in0);
	  REG_LOAD_VIR_32_64_COLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	  //load the second
	  DECLARE_REG_VIR_COLOR(in1);
	  REG_LOAD_VIR_32_64_COLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	  //merge the two and save
	  DECLARE_REG_VIR_COLOR(to_buf);
	  REG_VIR_COLOR_V0_MERGE(to_buf,in0,in1);
	  STORE_REG_VIR_32_64_COLOR(((VIR_32_64_COLOR*)send_buf)[idst_buf],to_buf);
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
	  DECLARE_REG_VIR_COLOR(in0);
	  REG_LOAD_VIR_32_64_COLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	  //load the second
	  DECLARE_REG_VIR_COLOR(in1);
	  REG_LOAD_VIR_32_64_COLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	  
	  //merge the two and save
	  DECLARE_REG_VIR_COLOR(to_buf);
	  REG_VIR_COLOR_V1_MERGE(to_buf,in0,in1);
	  STORE_REG_VIR_32_64_COLOR(((VIR_32_64_COLOR*)send_buf)[idst_buf],to_buf);
	}
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //pick data from vbuffer and put it in correct position
  THREADABLE_FUNCTION_0ARG(BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE)
  {
    GET_THREAD_ID();
    
    //short access
    const int fact=2;
    int v=vnode_paral_dir;
    VIR_32_64_COLOR *bgq_hopping_matrix_output_vdir_buffer=(VIR_32_64_COLOR*)send_buf+bord_volh/fact+8*loc_volh/fact;
    
    FORM_TWO_THREAD_TEAMS();
    
    //backward scatter (forward derivative)
    VIR_32_64_COLOR *base_out_bw=(VIR_32_64_COLOR*)send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
    VIR_32_64_COLOR *base_in_bw=bgq_hopping_matrix_output_vdir_buffer+0*vdir_bord_vol/fact;
    if(is_in_first_team)
      NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
	{
	  //load
	  DECLARE_REG_VIR_COLOR(reg_in);
	  REG_LOAD_VIR_32_64_COLOR(reg_in,base_in_bw[isrc]);
	  //transpose and store
	  DECLARE_REG_VIR_COLOR(reg_out);
	  REG_VIR_COLOR_TRANSPOSE(reg_out,reg_in);
	  STORE_REG_VIR_32_64_COLOR(base_out_bw[isrc*8+0+v],reg_out);
	}
    
    //forward scatter (backward derivative)
    VIR_32_64_COLOR *base_out_fw=(VIR_32_64_COLOR*)send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;  
    VIR_32_64_COLOR *base_in_fw=bgq_hopping_matrix_output_vdir_buffer+1*vdir_bord_vol/fact;
    if(is_in_second_team)
      NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
	{
	  //load
	  DECLARE_REG_VIR_COLOR(reg_in);
	  REG_LOAD_VIR_32_64_COLOR(reg_in,base_in_fw[isrc]);
	  //transpose and store
	  DECLARE_REG_VIR_COLOR(reg_out);
	  REG_VIR_COLOR_TRANSPOSE(reg_out,reg_in);
	  STORE_REG_VIR_32_64_COLOR(base_out_fw[isrc*8+4+v],reg_out);
	}
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //perform communications between VN and start all the communications between nodes
  THREADABLE_FUNCTION_0ARG(START_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS)
  {
    //shuffle data between virtual nodes and fill vdir out buffer
    if(paral_dir[vnode_paral_dir]) BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL();
    
    //start communications of scattered data to other nodes
    comm_start(EO_32_64_COLOR_COMM);
    
    //if v dir is not parallelized we have only to transpose between VN, and no real communication happens
    if(!paral_dir[vnode_paral_dir]) BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE();
  }
  THREADABLE_FUNCTION_END

  //finish the communications and put in place the communicated data
  THREADABLE_FUNCTION_1ARG(FINISH_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS, int,oe_or_eo)
  {
    GET_THREAD_ID();
    FORM_TWO_THREAD_TEAMS();
    
    //short access
    int v=vnode_paral_dir;
    const int fact=2;
    int *def_pos=viroe_or_vireo_hopping_matrix_output_pos[oe_or_eo].inter_fr_recv_pos;
    
    //wait communications end
    comm_wait(EO_32_64_COLOR_COMM);
    
    //vdir bw border (bw derivative): data goes to VN 0
    if(is_in_first_team)
      {
	if(paral_dir[v])
	  {
	    //inside incoming borders vdir is ordered naturally, while in the output data it comes first
	    VIR_32_64_COLOR *base_out=(VIR_32_64_COLOR*)send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;
	    VIR_32_64_COLOR *base_vdir_in=(VIR_32_64_COLOR*)send_buf+(bord_volh+8*loc_volh+1*vbord_vol/2)/fact;
	    VIR_32_64_COLOR *base_bord_in=(VIR_32_64_COLOR*)recv_buf+bord_offset[v]/2/fact;
	    
	    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/2/fact,thread_in_team_id,nthreads_in_team)
	      {
		//VN=0 must be filled with border
		DECLARE_REG_VIR_COLOR(in0);
		REG_LOAD_VIR_32_64_COLOR(in0,base_bord_in[isrc]);
		//VN=1 with buf0
		DECLARE_REG_VIR_COLOR(in1);
		REG_LOAD_VIR_32_64_COLOR(in1,base_vdir_in[2*isrc]);
		//merge and save
		DECLARE_REG_VIR_COLOR(to_dest);
		REG_VIR_COLOR_V0_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_COLOR(base_out[2*isrc*8+4+v],to_dest);
		
		//VN=1 with buf1
		REG_LOAD_VIR_32_64_COLOR(in1,base_vdir_in[2*isrc+1]);
		//merge and save
		REG_VIR_COLOR_V10_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_COLOR(base_out[(2*isrc+1)*8+4+v],to_dest);
	      }
	  }
	
	//other 3 bw borders
	for(int imu=0;imu<3;imu++)
	  {
	    int mu=perp_dir[v][imu];
	    VIR_32_64_COLOR *base_out=(VIR_32_64_COLOR*)send_buf+bord_volh/fact;
	    VIR_32_64_COLOR *base_in=(VIR_32_64_COLOR*)recv_buf;
	    NISSA_CHUNK_LOOP(isrc,bord_offset[mu]/2/fact,(bord_offset[mu]+bord_dir_vol[mu])/2/fact,
			     thread_in_team_id,nthreads_in_team)
	      SITE_COPY(base_out[def_pos[isrc]],base_in[isrc]);
	  }
      }
    
    //v fw border (fw derivative): data goes to VN 1
    if(is_in_second_team)
      {    
	if(paral_dir[v])
	  {
	    //inside incoming borders vdir is ordered naturally, while in the output data it comes first
	    VIR_32_64_COLOR *base_out=(VIR_32_64_COLOR*)send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
	    VIR_32_64_COLOR *base_vdir_in=(VIR_32_64_COLOR*)send_buf+(bord_volh+8*loc_volh+0*vbord_vol/2)/fact;
	    VIR_32_64_COLOR *base_bord_in=(VIR_32_64_COLOR*)recv_buf+(bord_vol/4+bord_offset[v]/2)/fact;
	    
	    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/2/fact,thread_in_team_id,nthreads_in_team)
	      {
		//VN=0 with buf1
		DECLARE_REG_VIR_COLOR(in0);
		REG_LOAD_VIR_32_64_COLOR(in0,base_vdir_in[2*isrc]);
		//VN=1 must be filled with border 0
		DECLARE_REG_VIR_COLOR(in1);
		REG_LOAD_VIR_32_64_COLOR(in1,base_bord_in[isrc]);
		//merge and save
		DECLARE_REG_VIR_COLOR(to_dest);
		REG_VIR_COLOR_V10_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_COLOR(base_out[2*isrc*8+0+v],to_dest);
		
		//VN=0 with buf1
		REG_LOAD_VIR_32_64_COLOR(in0,base_vdir_in[2*isrc+1]);
		//merge and save
		REG_VIR_COLOR_V1_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_COLOR(base_out[(2*isrc+1)*8+0+v],to_dest);
	      }
	  }
	
	//other 3 fw borders
	for(int imu=0;imu<3;imu++)
	  {
	    int mu=perp_dir[v][imu];
	    VIR_32_64_COLOR *base_out=(VIR_32_64_COLOR*)send_buf+bord_volh/fact;
	    VIR_32_64_COLOR *base_in=(VIR_32_64_COLOR*)recv_buf;
	    NISSA_CHUNK_LOOP(isrc,(bord_volh+bord_offset[mu])/2/fact,(bord_volh+bord_offset[mu]+bord_dir_vol[mu])/2/fact,
			     thread_in_team_id,nthreads_in_team)
	      SITE_COPY(base_out[def_pos[isrc]],base_in[isrc]);
	  }
      }
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
}

#undef PREC_TYPE
#undef VIR_32_64_COMPLEX
#undef VIR_32_64_COLOR
#undef VIR_32_64_SU3
#undef VIR_32_64_OCT_SU3
#undef EO_32_64_COLOR_COMM
#undef BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL
#undef BGQ_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE
#undef START_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS
#undef FINISH_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS
#undef HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D
#undef HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_STAGGERED_D_SUBTRACT_FROM_MASS2_TIMES_IN
#undef APPLY_32_64_STAGGERED_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM
#undef VIR_32_64_SU3_PREFETCH_NEXT
#undef REG_LOAD_VIR_32_64_COLOR_ADVANCING
#undef REG_LOAD_VIR_32_64_COLOR
#undef STORE_REG_VIR_32_64_COLOR
#undef STORE_REG_VIR_32_64_COLOR_ADVANCING
#undef VIR_32_64_COLOR_PREFETCH_NEXT
#undef REG_VIR_32_64_SU3_PROD_VIR_32_64_COLOR_LOAD_STORE
#undef REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_COLOR_LOAD_STORE
