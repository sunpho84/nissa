#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "bgq/bgq_macros.hpp"
#include "bgq/Wilson_hopping_matrix_lx_bgq.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ_bgq.hpp"
#include "new_types/complex.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

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
  
  THREADABLE_FUNCTION_5ARG(hopping_matrix_lx_tmclovQ_diag_term_bgq, bi_spincolor*,out, double,kappa, double,mu, bi_opt_as2t_su3*, Cl, bi_spincolor*,in)
  {
    GET_THREAD_ID();

    double A=-1/kappa,B=-2*mu;
    bi_complex diag[2]={{{+A,B},{+A,B}},{{-A,B},{-A,B}}};

    DECLARE_REG_BI_COMPLEX(reg_diag);

    //wait that all the terms are put in place
    THREAD_BARRIER();

    /* factor csw/2 already included in Cl
       +A  -B   0    0
       +B  -A   0    0
       0    0   +C  -D   BUT C and D get an additional
       0    0   +D  -C   sign -1 coming from g5
    */
	
    NISSA_PARALLEL_LOOP(i,0,loc_volh)
      {
	DECLARE_REG_BI_HALFSPINCOLOR(reg_out);
	DECLARE_REG_BI_SU3(U);
	DECLARE_REG_BI_HALFSPINCOLOR(reg_in);

	//
	
	//load first half of in
	REG_LOAD_BI_HALFSPINCOLOR(reg_in,in[i][0]);

	//multiply the first two terms
	REG_LOAD_BI_COMPLEX(reg_diag,diag[0]);
	REG_BI_HALFSPINCOLOR_PROD_COMPLEX(reg_out,reg_in,reg_diag);
	
	//deal with A
	REG_LOAD_BI_SU3(U,Cl[i][0]);
	REG_BI_SU3_SUMM_THE_PROD_BI_COLOR(reg_out_s0,U,reg_in_s0);
	REG_BI_SU3_SUBT_THE_PROD_BI_COLOR(reg_out_s1,U,reg_in_s1);
	REORDER_BARRIER();
	//deal with B
	REG_LOAD_BI_SU3(U,Cl[i][1]);
	REG_BI_SU3_SUBT_THE_PROD_BI_COLOR(reg_out_s0,U,reg_in_s1);
	REG_BI_SU3_SUMM_THE_PROD_BI_COLOR(reg_out_s1,U,reg_in_s0);
	
	STORE_REG_BI_HALFSPINCOLOR(out[i][0],reg_out);
	
	//

	//load second half of in
	REG_LOAD_BI_HALFSPINCOLOR(reg_in,in[i][2]);

	//multiply the other two terms
	REG_LOAD_BI_COMPLEX(reg_diag,diag[1]);
	REG_BI_HALFSPINCOLOR_PROD_COMPLEX(reg_out,reg_in,reg_diag);
	
	//deal with C
	REG_LOAD_BI_SU3(U,Cl[i][2]);
	REG_BI_SU3_SUBT_THE_PROD_BI_COLOR(reg_out_s0,U,reg_in_s0);
	REG_BI_SU3_SUMM_THE_PROD_BI_COLOR(reg_out_s1,U,reg_in_s1);
	REORDER_BARRIER();
	//deal with D
	REG_LOAD_BI_SU3(U,Cl[i][3]);
	REG_BI_SU3_SUMM_THE_PROD_BI_COLOR(reg_out_s0,U,reg_in_s1);
	REG_BI_SU3_SUBT_THE_PROD_BI_COLOR(reg_out_s1,U,reg_in_s0);
	
	STORE_REG_BI_HALFSPINCOLOR(out[i][2],reg_out);
      }

    //final sync
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_6ARG(apply_tmQ_bgq, bi_spincolor*,out, bi_oct_su3*,conf, double,kappa, double,mu, bi_opt_as2t_su3*,Cl, bi_spincolor*,in)
  {
    //compute on the surface and start communications
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(conf,0,vsurf_vol,in);
    start_Wilson_hopping_matrix_lx_bgq_communications();
    
    //compute on the bulk and put diagonal term while waiting to finish communications
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(conf,vsurf_vol,loc_volh,in);
    hopping_matrix_lx_tmclovQ_diag_term_bgq(out,kappa,mu,Cl,in);
    finish_Wilson_hopping_matrix_lx_bgq_communications();
    
    //put together all the 8 pieces
    hopping_matrix_lx_expand_to_Q_bgq(out);
  }
  THREADABLE_FUNCTION_END
}
