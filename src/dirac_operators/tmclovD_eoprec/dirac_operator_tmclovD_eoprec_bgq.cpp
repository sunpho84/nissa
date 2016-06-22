#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "bgq/bgq_macros.hpp"
#include "bgq/clover_term_bgq.hpp"
#include "bgq/intrinsic.hpp"
#include "bgq/Wilson_hopping_matrix_eo_or_oe_bgq.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.hpp"
#include "operations/su3_paths/clover_term.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  void inv_tmclovDee_or_oo_eos(vir_spincolor *out,vir_inv_clover_term_t *invCl,bool dag,vir_spincolor *in)
  {
    if(in==out) crash("in==out!");
    
    int high=0,low=1;
    
    GET_THREAD_ID();
    if(!dag)
      NISSA_PARALLEL_LOOP(X,0,loc_volh/2)
	{
	  unsafe_vir_halfspincolor_halfspincolor_times_vir_halfspincolor(&(out[X][2*high]),invCl[X][high],&(in[X][2*high]));
	  unsafe_vir_halfspincolor_halfspincolor_dag_times_vir_halfspincolor(&(out[X][2*low]),invCl[X][low],&(in[X][2*low]));
	}
    else
      NISSA_PARALLEL_LOOP(X,0,loc_volh/2)
	{
	  unsafe_vir_halfspincolor_halfspincolor_dag_times_vir_halfspincolor(&(out[X][2*high]),invCl[X][high],&(in[X][2*high]));
	  unsafe_vir_halfspincolor_halfspincolor_times_vir_halfspincolor(&(out[X][2*low]),invCl[X][low],&(in[X][2*low]));
	}
    
    set_borders_invalid(out);
  }
  
  //implement ee or oo part of Dirac operator, equation(3) and finish
  void minus_one_quarter_subt_tmclovDee_or_oo_eos(vir_spincolor *out,double kappa,vir_clover_term_t *Cl,bool dag,double mass,vir_spincolor *in)
  {
    if(dag) mass=-mass;
    
    if(in==out) crash("in==out!");
    
    //mass and kapap
    double a=1/(2*kappa),b=mass;
    vir_complex z[2]={{a,b, a,b},{a,-b, a,-b}};
    
    //coefficient for out
    double coef=-0.25;
    DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_OUT);
    REG_SPLAT_VIR_COMPLEX(VIR_REG_COEF_OUT,coef);
    
    //gamma5
    double g5[2]={1,-1};
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh/2)
      for(int high_low=0;high_low<2;high_low++)
	{
	  DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_IN);
	  REG_LOAD_VIR_COMPLEX(VIR_REG_COEF_IN,z[high_low]);
	  
	  DECLARE_REG_VIR_HALFSPINCOLOR(VIR_REG_OUT);
	  REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_OUT,out[ivol][NDIRAC/2*high_low]);
	  DECLARE_REG_VIR_HALFSPINCOLOR(VIR_REG_IN);
	  REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_IN,in[ivol][NDIRAC/2*high_low]);
	  
	  //multiply OUT{[0],[1]} by-+1/4
	  REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_OUT,out[ivol][2*high_low]);
	  REG_VIR_HALFSPINCOLOR_PROD_4DOUBLE(VIR_REG_OUT,VIR_REG_OUT,VIR_REG_COEF_OUT);
	  
	  //summ Z*in
	  REG_LOAD_VIR_HALFSPINCOLOR(VIR_REG_IN,in[ivol][2*high_low]);
	  REG_VIR_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(VIR_REG_OUT,VIR_REG_IN,VIR_REG_COEF_IN);
	  
	  REORDER_BARRIER();
	  
	  //terms related to A
	  DECLARE_REG_VIR_SU3(VIR_REG_Cl_0);
	  REG_LOAD_VIR_SU3(VIR_REG_Cl_0,Cl[ivol][0+2*high_low]);
	  REG_VIR_SU3_SUMM_THE_PROD_VIR_COLOR(VIR_REG_OUT_s0,VIR_REG_Cl_0,VIR_REG_IN_s0);
	  REG_VIR_SU3_SUBT_THE_PROD_VIR_COLOR(VIR_REG_OUT_s1,VIR_REG_Cl_0,VIR_REG_IN_s1);
	  
	  REORDER_BARRIER();
	  
	  //terms related to B
	  DECLARE_REG_VIR_SU3(VIR_REG_Cl_1);
	  REG_LOAD_VIR_SU3(VIR_REG_Cl_1,Cl[ivol][1+2*high_low]);
	  REG_VIR_SU3_DAG_SUMM_THE_PROD_VIR_COLOR(VIR_REG_OUT_s0,VIR_REG_Cl_1,VIR_REG_IN_s1);
	  REG_VIR_SU3_SUMM_THE_PROD_VIR_COLOR(VIR_REG_OUT_s1,VIR_REG_Cl_1,VIR_REG_IN_s0);
	  
	  REORDER_BARRIER();
	  
	  //gamma5
	  DECLARE_REG_VIR_COMPLEX(VIR_REG_COEF_G5);
	  REG_SPLAT_VIR_COMPLEX(VIR_REG_COEF_G5,g5[high_low]);
	  REG_VIR_HALFSPINCOLOR_PROD_4DOUBLE(VIR_REG_OUT,VIR_REG_OUT,VIR_REG_COEF_G5);
	  
	  STORE_REG_VIR_HALFSPINCOLOR(out[ivol][2*high_low],VIR_REG_OUT);
      }
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7)
  THREADABLE_FUNCTION_9ARG(tmclovDkern_eoprec_eos_bgq, vir_spincolor*,out, vir_spincolor*,temp, vir_oct_su3**,conf, double,kappa, vir_clover_term_t*,Cl_odd, vir_inv_clover_term_t*,invCl_evn, bool,dag, double,mass, vir_spincolor*,in)
  {
    if(in==out) crash("cannot work with in==out");
    
    tmn2Deo_eos_bgq(out,conf,in);
    inv_tmclovDee_or_oo_eos(temp,invCl_evn,dag,out);
    tmn2Doe_eos_bgq(out,conf,temp);
    
    minus_one_quarter_subt_tmclovDee_or_oo_eos(out,kappa,Cl_odd,dag,mass,in);
  }
  THREADABLE_FUNCTION_END
}
