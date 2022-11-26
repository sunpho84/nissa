#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/bench.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/gauge/Wilson_force.hpp"
#include "hmc/gauge/Symanzik_force.hpp"
#include "hmc/backfield.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#include "threads/threads.hpp"

#include "hmc/theory_pars.hpp"
#include "gluonic_action.hpp"

//#define DEBUG

namespace nissa
{
  /// Finish the computation multiplying for the conf and taking TA
  void gluonic_force_finish_computation(LxField<quad_su3>& F,
					const LxField<quad_su3>& conf)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[ivol][mu],F[ivol][mu]);
	  unsafe_su3_traceless_anti_hermitian_part(F[ivol][mu],temp);
	}
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Compute the gauge force
  void compute_gluonic_force_lx_conf_do_not_finish(LxField<quad_su3>& F,
				     const LxField<quad_su3>& conf,
				     const theory_pars_t& physics)
  {
    switch(physics.gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:
	Wilson_force_lx_conf(F,conf,physics.beta);
	break;
      case TLSYM_GAUGE_ACTION:
	Symanzik_force_lx_conf(F,conf,physics.beta,C1_TLSYM);
	break;
      case IWASAKI_GAUGE_ACTION:
	Symanzik_force_lx_conf(F,conf,physics.beta,C1_IWASAKI);
	break;
      default: crash("Unknown action");
      }
  }
  
  /// Take also the TA
  void compute_gluonic_force_lx_conf(LxField<quad_su3>& F,
				     const LxField<quad_su3>& conf,
				     const theory_pars_t& physics)
  {
    START_TIMING(gluon_force_time,ngluon_force);
    
#ifdef DEBUG
    vector_reset(F);
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori;
    gluonic_action(&act_ori,conf,physics->gauge_action_name,physics->beta);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<NCOL*NCOL-1;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus;
	gluonic_action(&act_minus,conf,physics->gauge_action_name,physics->beta);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus;
	gluonic_action(&act_plus,conf,physics->gauge_action_name,physics->beta);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    vector_reset(F);
#endif
    
    compute_gluonic_force_lx_conf_do_not_finish(F,conf,physics);
    
    //finish the calculation
    gluonic_force_finish_computation(F,conf);
    
#ifdef DEBUG
    master_printf("checking pure gauge force\n");
    master_printf("an\n");
    su3_print(F[0][0]);
    master_printf("nu\n");
    su3_print(nu);
    master_printf("nu_plus\n");
    su3_print(nu_plus);
    master_printf("nu_minus\n");
    su3_print(nu_minus);
    //crash("anna");
#endif
    
    //print the intensity of the force
    verbosity_lv2_master_printf("  Gluonic force average norm: %lg\n",sqrt(F.norm2()/glbVol));
    
    STOP_TIMING(gluon_force_time);
  }
}
