#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_vir.hpp"
#include "new_types/complex.hpp"
#include "Wilson_hopping_matrix_lx_bgq.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "bgq_macros.hpp"

#if (((defined VERSION_64) && (defined VERSION_32)) || ((!defined VERSION_64) && (!defined VERSION_32)))
 #error Impossible!
#endif

#if ((defined VERSION_64) && (!defined VERSION_32))
 #define PREC_TYPE double
 #define VIR_32_64_COMPLEX vir_complex
 #define VIR_32_64_HALFSPINCOLOR vir_halfspincolor
 #define VIR_32_64_SPINCOLOR vir_spincolor
 #define VIR_32_64_SU3 vir_su3
 #define VIR_32_64_OCT_SU3 vir_oct_su3
 #define EO_32_64_HALFSPINCOLOR_COMM eo_halfspincolor_comm
 #define BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL \
         bgq_double_Wilson_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill
 #define BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE \
         bgq_double_Wilson_hopping_matrix_oe_or_eo_vdir_VN_local_transpose
 #define START_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         start_double_Wilson_hopping_matrix_oe_or_eo_bgq_communications
 #define FINISH_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         finish_double_Wilson_hopping_matrix_oe_or_eo_bgq_communications
 #define HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_WILSON_2D \
         hopping_matrix_oe_or_eo_expand_to_double_Wilson_2D_bgq
#define APPLY_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM	\
         apply_double_Wilson_hopping_matrix_oe_or_eo_bgq_nocomm
 #define VIR_32_64_SU3_PREFETCH_NEXT VIR_SU3_PREFETCH_NEXT
 #define REG_LOAD_VIR_32_64_SPINCOLOR_ADVANCING REG_LOAD_VIR_SPINCOLOR_ADVANCING
 #define REG_LOAD_VIR_32_64_SPINCOLOR REG_LOAD_VIR_SPINCOLOR
 #define REG_LOAD_VIR_32_64_HALFSPINCOLOR REG_LOAD_VIR_HALFSPINCOLOR
 #define STORE_REG_VIR_32_64_SPINCOLOR STORE_REG_VIR_SPINCOLOR
 #define STORE_REG_VIR_32_64_SPINCOLOR_ADVANCING STORE_REG_VIR_SPINCOLOR_ADVANCING
 #define STORE_REG_VIR_32_64_HALFSPINCOLOR STORE_REG_VIR_HALFSPINCOLOR
 #define VIR_32_64_HALFSPINCOLOR_PREFETCH_NEXT VIR_HALFSPINCOLOR_PREFETCH_NEXT
 #define REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE REG_VIR_SU3_PROD_VIR_HALFSPINCOLOR_LOAD_STORE
 #define REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE REG_VIR_SU3_DAG_PROD_VIR_HALFSPINCOLOR_LOAD_STORE
 #define DER_TMQ_EXP_BGQ_32_64_HEADER DER_TMQ_EXP_BGQ_HEADER
#endif

#if ((defined VERSION_32) && (!defined VERSION_64))
 #define PREC_TYPE float
 #define VIR_32_64_COMPLEX vir_single_complex
 #define VIR_32_64_HALFSPINCOLOR vir_single_halfspincolor
 #define VIR_32_64_SPINCOLOR vir_single_spincolor
 #define VIR_32_64_SU3 vir_single_su3
 #define VIR_32_64_OCT_SU3 vir_single_oct_su3
 #define EO_32_64_HALFSPINCOLOR_COMM eo_single_halfspincolor_comm
 #define BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL \
         bgq_single_Wilson_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill
 #define BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE \
         bgq_single_Wilson_hopping_matrix_oe_or_eo_vdir_VN_local_transpose
 #define START_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         start_single_Wilson_hopping_matrix_oe_or_eo_bgq_communications
 #define FINISH_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS \
         finish_single_Wilson_hopping_matrix_oe_or_eo_bgq_communications
 #define HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_WILSON_2D \
         hopping_matrix_oe_or_eo_expand_to_single_Wilson_2D_bgq
#define APPLY_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM	\
         apply_single_Wilson_hopping_matrix_oe_or_eo_bgq_nocomm
 #define VIR_32_64_SU3_PREFETCH_NEXT VIR_SINGLE_SU3_PREFETCH_NEXT
 #define REG_LOAD_VIR_32_64_SPINCOLOR_ADVANCING REG_LOAD_VIR_SINGLE_SPINCOLOR_ADVANCING
 #define REG_LOAD_VIR_32_64_SPINCOLOR REG_LOAD_VIR_SINGLE_SPINCOLOR
 #define REG_LOAD_VIR_32_64_HALFSPINCOLOR REG_LOAD_VIR_SINGLE_HALFSPINCOLOR
 #define STORE_REG_VIR_32_64_SPINCOLOR STORE_REG_VIR_SINGLE_SPINCOLOR
 #define STORE_REG_VIR_32_64_SPINCOLOR_ADVANCING STORE_REG_VIR_SINGLE_SPINCOLOR_ADVANCING
 #define STORE_REG_VIR_32_64_HALFSPINCOLOR STORE_REG_VIR_SINGLE_HALFSPINCOLOR
 #define VIR_32_64_HALFSPINCOLOR_PREFETCH_NEXT VIR_SINGLE_HALFSPINCOLOR_PREFETCH_NEXT
 #define REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE REG_VIR_SINGLE_SU3_PROD_VIR_SINGLE_HALFSPINCOLOR_LOAD_STORE
 #define REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE REG_VIR_SINGLE_SU3_DAG_PROD_VIR_SINGLE_HALFSPINCOLOR_LOAD_STORE
 #define DER_TMQ_EXP_BGQ_32_64_HEADER DER_TMQ_EXP_BGQ_SINGLE_HEADER
#endif

#ifdef BGQ
 #define SITE_COPY(base_out,base_in)			\
  do							\
    {							\
      void *in=base_in,*out=base_out;			\
      DECLARE_REG_VIR_HALFSPINCOLOR(reg_temp);		\
      VIR_32_64_HALFSPINCOLOR_PREFETCH_NEXT(in);		\
      REG_LOAD_VIR_32_64_HALFSPINCOLOR(reg_temp,in);	\
      STORE_REG_VIR_32_64_HALFSPINCOLOR(out,reg_temp);	\
    }							\
  while(0)
#else
 #define SITE_COPY(out,in) VIR_32_64_HALFSPINCOLOR_COPY(out,in)
#endif

#define PROJ_HEADER(A)                          \
  do						\
    {						\
      REORDER_BARRIER();			\
      CACHE_PREFETCH(out+A);			\
      VIR_32_64_SU3_PREFETCH_NEXT(links[A]);	\
    }						\
  while(0)
  
namespace nissa
{
  //summ the eight contributions
  THREADABLE_FUNCTION_1ARG(HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_WILSON_2D, VIR_32_64_SPINCOLOR*,out)
  {
    GET_THREAD_ID();
    
    //result of split application
    VIR_32_64_HALFSPINCOLOR *bgq_hopping_matrix_output_data=(VIR_32_64_HALFSPINCOLOR*)send_buf+bord_volh/2;
    
    //define workload and point to the begin of each chunk
    NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,loc_volh/2,thread_id,NACTIVE_THREADS);
    
    //8 pieces
    for(int i=start;i<end;i++)
      {
	VIR_32_64_HALFSPINCOLOR *piece=bgq_hopping_matrix_output_data+i*8;
	
	DECLARE_REG_VIR_SPINCOLOR(reg_out);
	REG_SPLAT_VIR_SPINCOLOR(reg_out,0.0);
	
	DECLARE_REG_VIR_HALFSPINCOLOR(reg_temp);
	
	//TFW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[0]);
	REG_VIR_COLOR_SUMM(reg_out_s2,reg_out_s2,reg_temp_s0);
	REG_VIR_COLOR_SUMM(reg_out_s3,reg_out_s3,reg_temp_s1);
	
	//XFW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[1]);
	REG_VIR_COLOR_ISUBT(reg_out_s2,reg_out_s2,reg_temp_s1);
	REG_VIR_COLOR_ISUBT(reg_out_s3,reg_out_s3,reg_temp_s0);
	
	//YFW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[2]);
	REG_VIR_COLOR_SUBT(reg_out_s2,reg_out_s2,reg_temp_s1);
	REG_VIR_COLOR_SUMM(reg_out_s3,reg_out_s3,reg_temp_s0);
	
	//ZFW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[3]);
	REG_VIR_COLOR_ISUBT(reg_out_s2,reg_out_s2,reg_temp_s0);
	REG_VIR_COLOR_ISUMM(reg_out_s3,reg_out_s3,reg_temp_s1);
	
	//TBW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[4]);
	REG_VIR_COLOR_SUBT(reg_out_s2,reg_out_s2,reg_temp_s0);
	REG_VIR_COLOR_SUBT(reg_out_s3,reg_out_s3,reg_temp_s1);
	
	//XBW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[5]);
	REG_VIR_COLOR_ISUMM(reg_out_s2,reg_out_s2,reg_temp_s1);
	REG_VIR_COLOR_ISUMM(reg_out_s3,reg_out_s3,reg_temp_s0);
	
	//YBW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[6]);
	REG_VIR_COLOR_SUMM(reg_out_s2,reg_out_s2,reg_temp_s1);
	REG_VIR_COLOR_SUBT(reg_out_s3,reg_out_s3,reg_temp_s0);
	
	//ZBW
	DER_TMQ_EXP_BGQ_32_64_HEADER(reg_out,reg_temp,piece[7]);
	REG_VIR_COLOR_ISUMM(reg_out_s2,reg_out_s2,reg_temp_s0);
	REG_VIR_COLOR_ISUBT(reg_out_s3,reg_out_s3,reg_temp_s1);
	
	//store
        STORE_REG_VIR_32_64_SPINCOLOR(&(out[i]),reg_out);
      }
    
    //final sync
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_5ARG(APPLY_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM, VIR_32_64_OCT_SU3**,conf, int,istart, int,iend, VIR_32_64_SPINCOLOR*,in, int,oe_or_eo)
  {
    GET_THREAD_ID();
    
    VIR_32_64_HALFSPINCOLOR *out=(VIR_32_64_HALFSPINCOLOR*)send_buf;
    
    NISSA_PARALLEL_LOOP(i,istart,iend)
      {
	//take short access to link and output indexing
	int *iout=viroe_or_vireo_hopping_matrix_output_pos[oe_or_eo].inter_fr_in_pos+i*8;
	VIR_32_64_SU3 *links=(VIR_32_64_SU3*)(conf[oe_or_eo]+i);
	
	//declare
        DECLARE_REG_VIR_SPINCOLOR(reg_in);
        DECLARE_REG_VIR_HALFSPINCOLOR(reg_proj);
	
	//load
	REG_LOAD_VIR_32_64_SPINCOLOR(reg_in,in[i]);
	
        //T backward scatter (forward derivative)
        PROJ_HEADER(0);
        REG_VIR_COLOR_SUMM(reg_proj_s0,reg_in_s0,reg_in_s2);
        REG_VIR_COLOR_SUMM(reg_proj_s1,reg_in_s1,reg_in_s3);
        REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[0]],links[0],reg_proj);
        
        //X backward scatter (forward derivative)
        PROJ_HEADER(1);
        REG_VIR_COLOR_ISUMM(reg_proj_s0,reg_in_s0,reg_in_s3);
        REG_VIR_COLOR_ISUMM(reg_proj_s1,reg_in_s1,reg_in_s2);
        REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[1]],links[1],reg_proj);
        
        //Y backward scatter (forward derivative)
        PROJ_HEADER(2);
        REG_VIR_COLOR_SUMM(reg_proj_s0,reg_in_s0,reg_in_s3);
        REG_VIR_COLOR_SUBT(reg_proj_s1,reg_in_s1,reg_in_s2);
        REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[2]],links[2],reg_proj);
        
        //Z backward scatter (forward derivative)
        PROJ_HEADER(3);
        REG_VIR_COLOR_ISUMM(reg_proj_s0,reg_in_s0,reg_in_s2);
        REG_VIR_COLOR_ISUBT(reg_proj_s1,reg_in_s1,reg_in_s3);
        REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[3]],links[3],reg_proj);
        
        //T forward scatter (backward derivative)
        PROJ_HEADER(4);
        REG_VIR_COLOR_SUBT(reg_proj_s0,reg_in_s0,reg_in_s2);
        REG_VIR_COLOR_SUBT(reg_proj_s1,reg_in_s1,reg_in_s3);
        REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[4]],links[4],reg_proj);
        
        //X forward scatter (backward derivative)
        PROJ_HEADER(5);
        REG_VIR_COLOR_ISUBT(reg_proj_s0,reg_in_s0,reg_in_s3);
        REG_VIR_COLOR_ISUBT(reg_proj_s1,reg_in_s1,reg_in_s2);
        REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[5]],links[5],reg_proj);
        
        //Y forward scatter (backward derivative)
        PROJ_HEADER(6);
        REG_VIR_COLOR_SUBT(reg_proj_s0,reg_in_s0,reg_in_s3);
        REG_VIR_COLOR_SUMM(reg_proj_s1,reg_in_s1,reg_in_s2);
        REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[6]],links[6],reg_proj);
        
        //Z forward scatter (backward derivative)
        PROJ_HEADER(7);
        REG_VIR_COLOR_ISUBT(reg_proj_s0,reg_in_s0,reg_in_s2);
        REG_VIR_COLOR_ISUMM(reg_proj_s1,reg_in_s1,reg_in_s3);
        REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE(out[iout[7]],links[7],reg_proj);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //if virtual parallelized dir is really parallelized, fill send buffers
  THREADABLE_FUNCTION_0ARG(BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL)
  {
    GET_THREAD_ID();
    
    //form two thread team
    FORM_TWO_THREAD_TEAMS();
    
    //short access
    int v=vnode_paral_dir;
    const int fact=2;
    VIR_32_64_HALFSPINCOLOR *bgq_hopping_matrix_output_data=(VIR_32_64_HALFSPINCOLOR*)send_buf+bord_volh/fact; //half vol bisp = vol sp
    VIR_32_64_HALFSPINCOLOR *bgq_hopping_matrix_output_vdir_buffer=bgq_hopping_matrix_output_data+8*loc_volh/fact;
    
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
	  DECLARE_REG_VIR_HALFSPINCOLOR(in0);
	  REG_LOAD_VIR_32_64_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	  //load the second
	  DECLARE_REG_VIR_HALFSPINCOLOR(in1);
	  REG_LOAD_VIR_32_64_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	  //merge the two and save
	  DECLARE_REG_VIR_HALFSPINCOLOR(to_buf);
	  REG_VIR_HALFSPINCOLOR_V0_MERGE(to_buf,in0,in1);
	  STORE_REG_VIR_32_64_HALFSPINCOLOR(((VIR_32_64_HALFSPINCOLOR*)send_buf)[idst_buf],to_buf);
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
	  DECLARE_REG_VIR_HALFSPINCOLOR(in0);
	  REG_LOAD_VIR_32_64_HALFSPINCOLOR(in0,bgq_hopping_matrix_output_vdir_buffer[isrc]);
	  //load the second
	  DECLARE_REG_VIR_HALFSPINCOLOR(in1);
	  REG_LOAD_VIR_32_64_HALFSPINCOLOR(in1,bgq_hopping_matrix_output_vdir_buffer[isrc+1]);
	  
	  //merge the two and save
	  DECLARE_REG_VIR_HALFSPINCOLOR(to_buf);
	  REG_VIR_HALFSPINCOLOR_V1_MERGE(to_buf,in0,in1);
	  STORE_REG_VIR_32_64_HALFSPINCOLOR(((VIR_32_64_HALFSPINCOLOR*)send_buf)[idst_buf],to_buf);
	}
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //pick data from vbuffer and put it in correct position
  THREADABLE_FUNCTION_0ARG(BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE)
  {
    GET_THREAD_ID();
    
    //short access
    const int fact=2;
    int v=vnode_paral_dir;
    VIR_32_64_HALFSPINCOLOR *bgq_hopping_matrix_output_vdir_buffer=(VIR_32_64_HALFSPINCOLOR*)send_buf+bord_volh/fact+8*loc_volh/fact;
    
    FORM_TWO_THREAD_TEAMS();
    
    //backward scatter (forward derivative)
    VIR_32_64_HALFSPINCOLOR *base_out_bw=(VIR_32_64_HALFSPINCOLOR*)send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
    VIR_32_64_HALFSPINCOLOR *base_in_bw=bgq_hopping_matrix_output_vdir_buffer+0*vdir_bord_vol/fact;
    if(is_in_first_team)
      NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
	{
	  //load
	  DECLARE_REG_VIR_HALFSPINCOLOR(reg_in);
	  REG_LOAD_VIR_32_64_HALFSPINCOLOR(reg_in,base_in_bw[isrc]);
	  //transpose and store
	  DECLARE_REG_VIR_HALFSPINCOLOR(reg_out);
	  REG_VIR_HALFSPINCOLOR_TRANSPOSE(reg_out,reg_in);
	  STORE_REG_VIR_32_64_HALFSPINCOLOR(base_out_bw[isrc*8+0+v],reg_out);
	}
    
    //forward scatter (backward derivative)
    VIR_32_64_HALFSPINCOLOR *base_out_fw=(VIR_32_64_HALFSPINCOLOR*)send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;  
    VIR_32_64_HALFSPINCOLOR *base_in_fw=bgq_hopping_matrix_output_vdir_buffer+1*vdir_bord_vol/fact;
    if(is_in_second_team)
      NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/fact,thread_in_team_id,nthreads_in_team)
	{
	  //load
	  DECLARE_REG_VIR_HALFSPINCOLOR(reg_in);
	  REG_LOAD_VIR_32_64_HALFSPINCOLOR(reg_in,base_in_fw[isrc]);
	  //transpose and store
	  DECLARE_REG_VIR_HALFSPINCOLOR(reg_out);
	  REG_VIR_HALFSPINCOLOR_TRANSPOSE(reg_out,reg_in);
	  STORE_REG_VIR_32_64_HALFSPINCOLOR(base_out_fw[isrc*8+4+v],reg_out);
	}
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //perform communications between VN and start all the communications between nodes
  THREADABLE_FUNCTION_0ARG(START_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS)
  {
    //shuffle data between virtual nodes and fill vdir out buffer
    if(paral_dir[vnode_paral_dir]) BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL();
    
    //start communications of scattered data to other nodes
    comm_start(EO_32_64_HALFSPINCOLOR_COMM);
    
    //if v dir is not parallelized we have only to transpose between VN, and no real communication happens
    if(!paral_dir[vnode_paral_dir]) BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE();
  }
  THREADABLE_FUNCTION_END
  
  //finish the communications and put in place the communicated data
  THREADABLE_FUNCTION_1ARG(FINISH_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS, int,oe_or_eo)
  {
    GET_THREAD_ID();
    FORM_TWO_THREAD_TEAMS();
    
    //short access
    int v=vnode_paral_dir;
    const int fact=2;
    int *def_pos=viroe_or_vireo_hopping_matrix_output_pos[oe_or_eo].inter_fr_recv_pos;
    
    //wait communications end
    comm_wait(EO_32_64_HALFSPINCOLOR_COMM);
    
    //vdir bw border (bw derivative): data goes to VN 0
    if(is_in_first_team)
      {
	if(paral_dir[v])
	  {
	    //inside incoming borders vdir is ordered naturally, while in the output data it comes first
	    VIR_32_64_HALFSPINCOLOR *base_out=(VIR_32_64_HALFSPINCOLOR*)send_buf+(bord_volh+0*8*vdir_bord_vol)/fact;
	    VIR_32_64_HALFSPINCOLOR *base_vdir_in=(VIR_32_64_HALFSPINCOLOR*)send_buf+(bord_volh+8*loc_volh+1*vbord_vol/2)/fact;
	    VIR_32_64_HALFSPINCOLOR *base_bord_in=(VIR_32_64_HALFSPINCOLOR*)recv_buf+bord_offset[v]/2/fact;
	    
	    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/2/fact,thread_in_team_id,nthreads_in_team)
	      {
		//VN=0 must be filled with border
		DECLARE_REG_VIR_HALFSPINCOLOR(in0);
		REG_LOAD_VIR_32_64_HALFSPINCOLOR(in0,base_bord_in[isrc]);
		//VN=1 with buf0
		DECLARE_REG_VIR_HALFSPINCOLOR(in1);
		REG_LOAD_VIR_32_64_HALFSPINCOLOR(in1,base_vdir_in[2*isrc]);
		//merge and save
		DECLARE_REG_VIR_HALFSPINCOLOR(to_dest);
		REG_VIR_HALFSPINCOLOR_V0_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_HALFSPINCOLOR(base_out[2*isrc*8+4+v],to_dest);
		
		//VN=1 with buf1
		REG_LOAD_VIR_32_64_HALFSPINCOLOR(in1,base_vdir_in[2*isrc+1]);
		//merge and save
		REG_VIR_HALFSPINCOLOR_V10_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_HALFSPINCOLOR(base_out[(2*isrc+1)*8+4+v],to_dest);
	      }
	  }
	
	//other 3 bw borders
	for(int imu=0;imu<3;imu++)
	  {
	    int mu=perp_dir[v][imu];
	    VIR_32_64_HALFSPINCOLOR *base_out=(VIR_32_64_HALFSPINCOLOR*)send_buf+bord_volh/fact;
	    VIR_32_64_HALFSPINCOLOR *base_in=(VIR_32_64_HALFSPINCOLOR*)recv_buf;
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
	    VIR_32_64_HALFSPINCOLOR *base_out=(VIR_32_64_HALFSPINCOLOR*)send_buf+(bord_volh+1*8*vdir_bord_vol)/fact;
	    VIR_32_64_HALFSPINCOLOR *base_vdir_in=(VIR_32_64_HALFSPINCOLOR*)send_buf+(bord_volh+8*loc_volh+0*vbord_vol/2)/fact;
	    VIR_32_64_HALFSPINCOLOR *base_bord_in=(VIR_32_64_HALFSPINCOLOR*)recv_buf+(bord_vol/4+bord_offset[v]/2)/fact;
	    
	    NISSA_CHUNK_LOOP(isrc,0,vdir_bord_vol/2/fact,thread_in_team_id,nthreads_in_team)
	      {
		//VN=0 with buf1
		DECLARE_REG_VIR_HALFSPINCOLOR(in0);
		REG_LOAD_VIR_32_64_HALFSPINCOLOR(in0,base_vdir_in[2*isrc]);
		//VN=1 must be filled with border 0
		DECLARE_REG_VIR_HALFSPINCOLOR(in1);
		REG_LOAD_VIR_32_64_HALFSPINCOLOR(in1,base_bord_in[isrc]);
		//merge and save
		DECLARE_REG_VIR_HALFSPINCOLOR(to_dest);
		REG_VIR_HALFSPINCOLOR_V10_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_HALFSPINCOLOR(base_out[2*isrc*8+0+v],to_dest);
		
		//VN=0 with buf1
		REG_LOAD_VIR_32_64_HALFSPINCOLOR(in0,base_vdir_in[2*isrc+1]);
		//merge and save
		REG_VIR_HALFSPINCOLOR_V1_MERGE(to_dest,in0,in1);
		STORE_REG_VIR_32_64_HALFSPINCOLOR(base_out[(2*isrc+1)*8+0+v],to_dest);
	      }
	  }
	
	//other 3 fw borders
	for(int imu=0;imu<3;imu++)
	  {
	    int mu=perp_dir[v][imu];
	    VIR_32_64_HALFSPINCOLOR *base_out=(VIR_32_64_HALFSPINCOLOR*)send_buf+bord_volh/fact;
	    VIR_32_64_HALFSPINCOLOR *base_in=(VIR_32_64_HALFSPINCOLOR*)recv_buf;
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
#undef VIR_32_64_SPINCOLOR
#undef VIR_32_64_HALFSPINCOLOR
#undef VIR_32_64_SU3
#undef VIR_32_64_OCT_SU3
#undef EO_32_64_HALFSPINCOLOR_COMM
#undef BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_COMM_AND_BUFF_FILL
#undef BGQ_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_VDIR_VN_LOCAL_TRANSPOSE
#undef START_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS
#undef FINISH_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_COMMUNICATIONS
#undef HOPPING_MATRIX_OE_OR_EO_EXPAND_TO_32_64_WILSON_2D
#undef APPLY_32_64_WILSON_HOPPING_MATRIX_OE_OR_EO_BGQ_NOCOMM
#undef VIR_32_64_SU3_PREFETCH_NEXT
#undef REG_LOAD_VIR_32_64_SPINCOLOR_ADVANCING
#undef REG_LOAD_VIR_32_64_SPINCOLOR
#undef STORE_REG_VIR_32_64_SPINCOLOR
#undef STORE_REG_VIR_32_64_SPINCOLOR_ADVANCING
#undef REG_LOAD_VIR_32_64_HALFSPINCOLOR
#undef STORE_REG_VIR_32_64_HALFSPINCOLOR
#undef VIR_32_64_HALFSPINCOLOR_PREFETCH_NEXT
#undef REG_VIR_32_64_SU3_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE
#undef DER_TMQ_EXP_BGQ_32_64_HEADER
#undef REG_VIR_32_64_SU3_DAG_PROD_VIR_32_64_HALFSPINCOLOR_LOAD_STORE
