#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "bgq/bgq_macros.hpp"
#include "bgq/Wilson_hopping_matrix_lx_bgq.hpp"
#include "communicate/communicate.hpp"
#include "new_types/complex.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  /*
    application of hopping matrix goes among the following steps:
    1) hopping matrix on the surf
    2) start communications
    3) hopping matrix on the bulk
    4) finish communications
    then the application of Q requires to expand the halfspincolor into full spincolor and summ the diagonal part
  */
  
  THREADABLE_FUNCTION_4ARG(hopping_matrix_lx_tmQ_diag_term_bgq, vir_spincolor*,out, double,kappa, double,mu, vir_spincolor*,in)
  {
    GET_THREAD_ID();
    
    double A=-1/kappa,B=-2*mu;
    vir_complex diag[2]={{{+A,B},{+A,B}},{{-A,B},{-A,B}}};
    
#ifdef BGQ
    DECLARE_REG_VIR_COMPLEX(reg_diag_0);
    DECLARE_REG_VIR_COMPLEX(reg_diag_1);
    REG_LOAD_VIR_COMPLEX(reg_diag_0,diag[0]);
    REG_LOAD_VIR_COMPLEX(reg_diag_1,diag[1]);
#endif
    
    //wait that all the terms are put in place
    THREAD_BARRIER();
    
    NISSA_PARALLEL_LOOP(i,0,loc_volh)
      {
#ifdef BGQ
	DECLARE_REG_VIR_SPINCOLOR(reg_out);
	
	//multiply in by the diagonal part (apart from the final term -0.5 added separately when summing non-diag)
	DECLARE_REG_VIR_SPINCOLOR(reg_in);
	REG_LOAD_VIR_SPINCOLOR(reg_in,in[i]);
	REG_VIR_COLOR_PROD_COMPLEX(reg_out_s0,reg_in_s0,reg_diag_0);
	REG_VIR_COLOR_PROD_COMPLEX(reg_out_s1,reg_in_s1,reg_diag_0);
	REG_VIR_COLOR_PROD_COMPLEX(reg_out_s2,reg_in_s2,reg_diag_1);
	VIR_SPINCOLOR_PREFETCH_NEXT(in[i]);
	REG_VIR_COLOR_PROD_COMPLEX(reg_out_s3,reg_in_s3,reg_diag_1);
	
	STORE_REG_VIR_SPINCOLOR(&(out[i]),reg_out);
#else
	vir_spincolor temp;
	
	//multiply in by the diagonal part
	DIAG_TMQ(out[i],diag,in[i]);
#endif
      }
    NISSA_PARALLEL_LOOP_END;
    
    //final sync
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_1ARG(hopping_matrix_lx_expand_to_Q_bgq, vir_spincolor*,out)
  {
    GET_THREAD_ID();
    
    //wait that all the terms are put in place
    THREAD_BARRIER();

#ifdef BGQ
    vir_complex mone_half={{-0.5,-0.5},{-0.5,-0.5}};
    DECLARE_REG_VIR_COMPLEX(reg_mone_half);
    REG_LOAD_VIR_COMPLEX(reg_mone_half,mone_half);
#endif
    
    vir_halfspincolor *bgq_hopping_matrix_output_data=(vir_halfspincolor*)send_buf+bord_volh;
    
    NISSA_PARALLEL_LOOP(i,0,loc_volh)
      {
#ifdef BGQ
	DECLARE_REG_VIR_SPINCOLOR(reg_out);
	REG_LOAD_VIR_SPINCOLOR(reg_out,out[i]);
	
	//8 pieces
	{
	  vir_halfspincolor *piece=bgq_hopping_matrix_output_data+i*8;
	  DECLARE_REG_VIR_HALFSPINCOLOR(reg_temp);
	  
	  //TFW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[0]);
	  REG_VIR_COLOR_SUBT(reg_out_s2,reg_out_s2,reg_temp_s0);
	  REG_VIR_COLOR_SUBT(reg_out_s3,reg_out_s3,reg_temp_s1);
	  
	  //XFW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[1]);
	  REG_VIR_COLOR_ISUMM(reg_out_s2,reg_out_s2,reg_temp_s1);
	  REG_VIR_COLOR_ISUMM(reg_out_s3,reg_out_s3,reg_temp_s0);
	  
	  //YFW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[2]);
	  REG_VIR_COLOR_SUMM(reg_out_s2,reg_out_s2,reg_temp_s1);
	  REG_VIR_COLOR_SUBT(reg_out_s3,reg_out_s3,reg_temp_s0);
	  
	  //ZFW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[3]);
	  REG_VIR_COLOR_ISUMM(reg_out_s2,reg_out_s2,reg_temp_s0);
	  REG_VIR_COLOR_ISUBT(reg_out_s3,reg_out_s3,reg_temp_s1);
	  
	  //TBW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[4]);
	  REG_VIR_COLOR_SUMM(reg_out_s2,reg_out_s2,reg_temp_s0);
	  REG_VIR_COLOR_SUMM(reg_out_s3,reg_out_s3,reg_temp_s1);
	  
	  //XBW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[5]);
	  REG_VIR_COLOR_ISUBT(reg_out_s2,reg_out_s2,reg_temp_s1);
	  REG_VIR_COLOR_ISUBT(reg_out_s3,reg_out_s3,reg_temp_s0);
	  
	  //YBW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[6]);
	  REG_VIR_COLOR_SUBT(reg_out_s2,reg_out_s2,reg_temp_s1);
	  REG_VIR_COLOR_SUMM(reg_out_s3,reg_out_s3,reg_temp_s0);
	  
	  VIR_SPINCOLOR_PREFETCH_NEXT(out[i]);
	  
	  //ZBW
	  DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_temp,piece[7]);
	  REG_VIR_COLOR_ISUBT(reg_out_s2,reg_out_s2,reg_temp_s0);
	  REG_VIR_COLOR_ISUMM(reg_out_s3,reg_out_s3,reg_temp_s1);
	}
	
	//put final -0.5
	REG_VIR_SPINCOLOR_PROD_4DOUBLE(reg_out,reg_out,reg_mone_half);
	STORE_REG_VIR_SPINCOLOR(&(out[i]),reg_out);
#else
	//summ 8 contributions
	TFW_DER_TMQ_EXP(out[i],piece[0]);
	XFW_DER_TMQ_EXP(out[i],piece[1]);
	YFW_DER_TMQ_EXP(out[i],piece[2]);
	ZFW_DER_TMQ_EXP(out[i],piece[3]);
	TBW_DER_TMQ_EXP(out[i],piece[4]);
	XBW_DER_TMQ_EXP(out[i],piece[5]);
	YBW_DER_TMQ_EXP(out[i],piece[6]);
	ZBW_DER_TMQ_EXP(out[i],piece[7]);
	
	//put final -0.5
	VIR_SPINCOLOR_PROD_DOUBLE(out[i],out[i],-0.5);
#endif
      }
    NISSA_PARALLEL_LOOP_END;
    
    //final sync
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_5ARG(apply_tmQ_bgq, vir_spincolor*,out, vir_oct_su3*,conf, double,kappa, double,mu, vir_spincolor*,in)
  {
    //compute on the surface and start communications
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(conf,0,vsurf_vol,in);
    start_Wilson_hopping_matrix_lx_bgq_communications();
    
    //compute on the bulk and multiply the diagonal part while waiting to finish communications
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(conf,vsurf_vol,loc_volh,in);
    hopping_matrix_lx_tmQ_diag_term_bgq(out,kappa,mu,in);
    finish_Wilson_hopping_matrix_lx_bgq_communications();
    
    //put together all the 8 pieces
    hopping_matrix_lx_expand_to_Q_bgq(out);
  }
  THREADABLE_FUNCTION_END
}
