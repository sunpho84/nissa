#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/theory_pars.hpp"
#include "new_types/su3.hpp"
#include "measures/gauge/topological_charge.hpp"
#include "operations/smearing/stout.hpp"
#include "threads/threads.hpp"
#include "routines/ios.hpp"

#include "gluonic_force.hpp"

//#define DEBUG

#ifdef DEBUG
 #include "topological_action.hpp"
#endif

namespace nissa
{
  //compute the topodynamical potential
  double compute_topodynamical_potential_der(const topotential_pars_t *pars,
					     const LxField<quad_su3>& conf)
  {
    crash("reimplement");//topological_charge.cpp not compiled
    // const double Q=total_topological_charge_lx_conf(conf);
    
    // return pars->compute_pot_der(Q);
    return {};
  }
  
  //common part, for staples and potential if needed
  void compute_topological_force_lx_conf_internal(LxField<quad_su3>& F,
						  const LxField<quad_su3>& conf,
						  const topotential_pars_t *pars)
  {
    crash("reimplement");//topological_charge.cpp not compiled
    //compute the staples
    // topological_staples(F,conf);
    
    // //compute the potential
    // double pot=0;
    // switch(pars->flag)
    //   {
    //   case 1:
    // 	pot=pars->theta;
    // 	break;
    //   case 2:
    // 	pot=compute_topodynamical_potential_der(pars,conf);
    // 	break;
    //   default:
    // 	crash("unknown way to compute topological potential %d",pars->flag);
    //   }
    
    // //normalize
    // const double norm=pot/(M_PI*M_PI*128);
    // F*=norm;
  }
  
  //compute the topological force
  void compute_topological_force_lx_conf(LxField<quad_su3>& F,
					 LxField<quad_su3>& conf,
					 const topotential_pars_t* pars)
  {
    crash("reimplement");
//     verbosity_lv1_master_printf("Computing topological force\n");
    
//     //compute the staples
//     if(pars->stout_pars.nlevels==0) compute_topological_force_lx_conf_internal(F,conf,pars);
//     else
//       {
// 	//allocate the stack of confs: conf is binded to sme_conf[0]
// 	std::vector<LxField<quad_su3>*> sme_conf;
//         stout_smear_conf_stack_allocate(sme_conf,conf,pars->stout_pars.nlevels);
        
//         //smear iteratively retaining all the stack
//         stout_smear_whole_stack(sme_conf,conf,&(pars->stout_pars));
        
//         //compute the force in terms of the most smeared conf
// 	compute_topological_force_lx_conf_internal(F,*sme_conf[pars->stout_pars.nlevels],pars);
	
//         //remap the force backward
//         stouted_force_remap(F,sme_conf,&(pars->stout_pars));
	
// 	//now free the stack of confs
//         stout_smear_conf_stack_free(sme_conf);
//       }
    
//     //take TA
//     gluonic_force_finish_computation(F,conf);
    
// #ifdef DEBUG
//     double eps=1e-5;
    
//     //store initial link and compute action
//     su3 sto;
//     su3_copy(sto,conf[0][0]);
//     double act_ori=topotential_action(conf,*pars);
    
//     //store derivative
//     su3 nu_plus,nu_minus;
//     su3_put_to_zero(nu_plus);
//     su3_put_to_zero(nu_minus);
    
//     for(int igen=0;igen<NCOL*NCOL-1;igen++)
//       {
// 	//prepare increment and change
// 	su3 ba;
// 	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	
// 	su3 exp_mod;
// 	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
// 	//change -, compute action
// 	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
// 	double act_minus=topotential_action(conf,*pars);
	
// 	//change +, compute action
// 	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
// 	double act_plus=topotential_action(conf,*pars);
	
// 	//set back everything
// 	su3_copy(conf[0][0],sto);
	
// 	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
// 	double gr_plus=-(act_plus-act_ori)/eps;
// 	double gr_minus=-(act_ori-act_minus)/eps;
// 	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
// 	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
//       }
    
//     //take the average
//     su3 nu;
//     su3_summ(nu,nu_plus,nu_minus);
//     su3_prodassign_double(nu,0.5);
    
//     master_printf("Comparing topotential force\n");
//     master_printf("an\n");
//     su3_print(F[0][0]);
//     master_printf("nu\n");
//     su3_print(nu);
//     master_printf("nu_plus\n");
//     su3_print(nu_plus);
//     master_printf("nu_minus\n");
//     su3_print(nu_minus);
//     crash("anna");
// #endif
  }
}
