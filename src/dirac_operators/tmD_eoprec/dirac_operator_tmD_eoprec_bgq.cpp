#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "bgq/bgq_macros.hpp"
#include "bgq/intrinsic.hpp"
#include "threads/threads.hpp"

#include "dirac_operator_tmD_eoprec_bgq.hpp"

namespace nissa
{
  void inv_tmDee_or_oo_eos(vir_spincolor *out,double kappa,double mass,vir_spincolor *in)
  {
    double a=1/(2*kappa),b=mass,nrm=a*a+b*b;
    a/=nrm;
    b/=-nrm;
    vir_complex z={a,b,a,b};
    DECLARE_REG_VIR_COMPLEX(VIR_REG_Z);
    REG_LOAD_VIR_COMPLEX(VIR_REG_Z,z);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh/2)
      {
	DECLARE_REG_VIR_HALFSPINCOLOR(VIR_REG_IN);
	DECLARE_REG_VIR_HALFSPINCOLOR(VIR_REG_OUT);
	
	REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_IN,in[X][0]);
	REG_VIR_HALFSPINCOLOR_PROD_COMPLEX(VIR_REG_OUT,VIR_REG_IN,VIR_REG_Z);
	STORE_REG_VIR_HALFSPINCOLOR(out[X][0],VIR_REG_OUT);
	
	REORDER_BARRIER();
	
	REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_IN,in[X][2]);
	REG_VIR_HALFSPINCOLOR_PROD_CONJ2_COMPLEX(VIR_REG_OUT,VIR_REG_IN,VIR_REG_Z);
	STORE_REG_VIR_HALFSPINCOLOR(out[X][2],VIR_REG_OUT);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement ee or oo part of Dirac operator, equation(3) and finish
  void minus_one_quarter_subt_tmD_or_Qee_or_oo_eos(vir_spincolor *out,double kappa,double mass,vir_spincolor *in,int D_or_Q)
  {
    DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_OUT_01);
    REG_SPLAT_VIR_COMPLEX(VIR_REG_COEF_OUT_01,-0.25);
    DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_OUT_23); //-+1 according to D or Q
    if(D_or_Q==0) REG_SPLAT_VIR_COMPLEX(VIR_REG_COEF_OUT_23,+0.25);
    else          REG_SPLAT_VIR_COMPLEX(VIR_REG_COEF_OUT_23,-0.25);
    
    double a=1/(2*kappa),b=mass;
    vir_complex z={a,b, a,b};
    DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_IN_01);
    REG_LOAD_VIR_COMPLEX(VIR_REG_COEF_IN_01,z);
    vir_complex zc={a,-b, a,-b}; //complex conjugate of Z (when using Q)
    vir_complex mzc={-a,b, -a,b}; //minus complex conjugate of Z (when using D)
    DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_IN_23);
    if(D_or_Q==0) REG_LOAD_VIR_COMPLEX(VIR_REG_COEF_IN_23,mzc);
    else          REG_LOAD_VIR_COMPLEX(VIR_REG_COEF_IN_23,zc);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh/2)
      {
	DECLARE_REG_VIR_HALFSPINCOLOR(VIR_REG_OUT);
	DECLARE_REG_VIR_HALFSPINCOLOR(VIR_REG_IN);
	
	//multiply OUT{[0],[1]} by -1/4
	REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_OUT,out[ivol][0]);
	REG_VIR_HALFSPINCOLOR_PROD_4DOUBLE(VIR_REG_OUT,VIR_REG_OUT,VIR_REG_COEF_OUT_01);
	
	//summ Z*in
	REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_IN,in[ivol][0]);
	REG_VIR_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(VIR_REG_OUT,VIR_REG_IN,VIR_REG_COEF_IN_01);
        STORE_REG_VIR_HALFSPINCOLOR(out[ivol][0],VIR_REG_OUT);
	
	REORDER_BARRIER();
	
	//multiply OUT{[2],[3]} by +1/4
	REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_OUT,out[ivol][2]);
	REG_VIR_HALFSPINCOLOR_PROD_4DOUBLE(VIR_REG_OUT,VIR_REG_OUT,VIR_REG_COEF_OUT_23);
	
	//subt (Z^+)*in
	REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_IN,in[ivol][2]);
	REG_VIR_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(VIR_REG_OUT,VIR_REG_IN,VIR_REG_COEF_IN_23);
        STORE_REG_VIR_HALFSPINCOLOR(out[ivol][2],VIR_REG_OUT);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7)
  THREADABLE_FUNCTION_6ARG(tmD_or_Qkern_eoprec_eos_bgq, vir_spincolor*,out, vir_oct_su3**,conf, double,kappa, double,mass, vir_spincolor*,in, int,D_or_Q)
  {
    if(in==out) crash("cannot work with in==out");
    
    tmn2Deo_eos_bgq(out,conf,in);
    inv_tmDee_or_oo_eos(out,kappa,mass,out);
    tmn2Doe_eos_bgq(out,conf,out);
    
    minus_one_quarter_subt_tmD_or_Qee_or_oo_eos(out,kappa,mass,in,D_or_Q);
  }
  THREADABLE_FUNCTION_END
}
