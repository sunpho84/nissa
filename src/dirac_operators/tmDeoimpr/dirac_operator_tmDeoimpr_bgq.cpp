#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/macros.hpp"
#include "base/thread_macros.hpp"
#include "bgq/bgq_macros.hpp"
#include "bgq/intrinsic.hpp"
#include "bgq/Wilson_hopping_matrix_eo_or_oe_bgq.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  void inv_tmDee_or_oo_eos(bi_spincolor *out,double kappa,double mass,bi_spincolor *in)
  {
    double a=1/(2*kappa),b=mass,nrm=a*a+b*b;
    a/=nrm;
    b/=-nrm;
    bi_complex z={a,b,a,b};
    DECLARE_REG_BI_COMPLEX(BI_REG_Z);
    REG_LOAD_BI_COMPLEX(BI_REG_Z,z);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh/2)
      {
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_IN);
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_OUT);
	
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[X][0]);
	REG_BI_HALFSPINCOLOR_PROD_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_Z);
	STORE_REG_BI_HALFSPINCOLOR(out[X][0],BI_REG_OUT);
	
	REORDER_BARRIER();
	
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[X][2]);
	REG_BI_HALFSPINCOLOR_PROD_CONJ2_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_Z);
	STORE_REG_BI_HALFSPINCOLOR(out[X][2],BI_REG_OUT);
      }
    
    set_borders_invalid(out);
  }
  
  //implement ee or oo part of Dirac operator, equation(3) and finish
  void minus_one_quarter_subt_tmD_or_Qee_or_oo_eos(bi_spincolor *out,double kappa,double mass,bi_spincolor *in,int D_or_Q)
  {
    DECLARE_REG_BI_COMPLEX(BI_REG_COEF_OUT_01);
    REG_SPLAT_BI_COMPLEX(BI_REG_COEF_OUT_01,-0.25);
    DECLARE_REG_BI_COMPLEX(BI_REG_COEF_OUT_23); //-+1 according to D or Q
    if(D_or_Q==0) REG_SPLAT_BI_COMPLEX(BI_REG_COEF_OUT_23,+0.25);
    else          REG_SPLAT_BI_COMPLEX(BI_REG_COEF_OUT_23,-0.25);
    
    double a=1/(2*kappa),b=mass;
    bi_complex z={a,b, a,b};
    DECLARE_REG_BI_COMPLEX(BI_REG_COEF_IN_01);
    REG_LOAD_BI_COMPLEX(BI_REG_COEF_IN_01,z);
    bi_complex zc={a,-b, a,-b}; //complex conjugate of Z (when using Q)
    bi_complex mzc={-a,b, -a,b}; //minus complex conjugate of Z (when using D)
    DECLARE_REG_BI_COMPLEX(BI_REG_COEF_IN_23);
    if(D_or_Q==0) REG_LOAD_BI_COMPLEX(BI_REG_COEF_IN_23,mzc);
    else          REG_LOAD_BI_COMPLEX(BI_REG_COEF_IN_23,zc);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh/2)
      {
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_OUT);
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_IN);
	
	//multiply OUT{[0],[1]} by -1/4
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_OUT,out[ivol][0]);
	REG_BI_HALFSPINCOLOR_PROD_4DOUBLE(BI_REG_OUT,BI_REG_OUT,BI_REG_COEF_OUT_01);
	
	//summ Z*in
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[ivol][0]);
	REG_BI_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_COEF_IN_01);
        STORE_REG_BI_HALFSPINCOLOR(out[ivol][0],BI_REG_OUT);
	
	REORDER_BARRIER();
	
	//multiply OUT{[2],[3]} by +1/4
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_OUT,out[ivol][2]);
	REG_BI_HALFSPINCOLOR_PROD_4DOUBLE(BI_REG_OUT,BI_REG_OUT,BI_REG_COEF_OUT_23);
	
	//subt (Z^+)*in
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[ivol][2]);
	REG_BI_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_COEF_IN_23);
        STORE_REG_BI_HALFSPINCOLOR(out[ivol][2],BI_REG_OUT);
      }
    set_borders_invalid(out);
  }
  
  void tmn2Doe_or_tmn2Deo_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,int oe_or_eo,bi_spincolor *in)
  {
    //compute on the surface and start communications
    apply_double_Wilson_hopping_matrix_oe_or_eo_bgq_nocomm(conf,0,vsurf_volh,in,oe_or_eo);
    start_double_Wilson_hopping_matrix_oe_or_eo_bgq_communications();
    
    //compute on the bulk and finish communications
    apply_double_Wilson_hopping_matrix_oe_or_eo_bgq_nocomm(conf,vsurf_volh,loc_volh/2,in,oe_or_eo);
    finish_double_Wilson_hopping_matrix_oe_or_eo_bgq_communications(oe_or_eo);
    
    //put together all the 8 pieces
    hopping_matrix_oe_or_eo_expand_to_double_Wilson_2D_bgq(out);
  }
  
  //wrappers
  void tmn2Doe_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,bi_spincolor *in){tmn2Doe_or_tmn2Deo_eos_bgq(out,conf,0,in);}
  void tmn2Deo_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,bi_spincolor *in){tmn2Doe_or_tmn2Deo_eos_bgq(out,conf,1,in);}
  
  //implement Koo defined in equation (7)
  THREADABLE_FUNCTION_6ARG(tmD_or_Qkern_eoprec_eos_bgq, bi_spincolor*,out, bi_oct_su3**,conf, double,kappa, double,mass, bi_spincolor*,in, int,D_or_Q)
  {
    if(in==out) crash("cannot work with in==out");
    
    tmn2Deo_eos_bgq(out,conf,in);
    inv_tmDee_or_oo_eos(out,kappa,mass,out);
    tmn2Doe_eos_bgq(out,conf,out);
    
    minus_one_quarter_subt_tmD_or_Qee_or_oo_eos(out,kappa,mass,in,D_or_Q);
  }
  THREADABLE_FUNCTION_END
}
