#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_lx.hpp"
 #include "geometry/geometry_mix.hpp"
 #include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "operations/smearing/stout.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#include "routines/ios.hpp"

#include "gluonic_force.hpp"

namespace nissa
{
  //compute the topodynamical potential
  double compute_topodynamical_potential_der(topotential_pars_t *pars,quad_su3 *conf)
  {
    double Q;
    total_topological_charge_lx_conf(&Q,conf);
    
    return pars->compute_pot_der(Q);
  }
  
  //common part, for staples and potential if needed
  void compute_topological_force_lx_conf_internal(quad_su3 *F,quad_su3 *conf,topotential_pars_t *pars,bool phase_pres)
  {
    GET_THREAD_ID();
    
    //compute the staples
    topological_staples(F,conf);
    
    //compute the potential
    double pot=0;
    switch(pars->flag)
      {
      case 1: pot=pars->theta;break;
      case 2: pot=compute_topodynamical_potential_der(pars,conf); break;
      default: crash("unknown way to compute topological potential %d",pars->flag);
      }
    
    //normalize
    double norm=pot/(M_PI*M_PI*128);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	safe_su3_hermitian_prod_double(F[ivol][mu],F[ivol][mu],norm);
    set_borders_invalid(F);
    
    //add the stag phases to the force term, to cancel the one entering the force
    if(phase_pres) addrem_stagphases_to_lx_conf(F);
  }
  
  //compute the topological force
  THREADABLE_FUNCTION_4ARG(compute_topological_force_lx_conf, quad_su3*,F, quad_su3*,conf, topotential_pars_t*,pars, bool,phase_pres)
  {
    verbosity_lv1_master_printf("Computing topological force\n");
    
    //compute the staples
    if(pars->stout_pars.nlevls==0) compute_topological_force_lx_conf_internal(F,conf,pars,phase_pres);
    else
      {
	//allocate the stack of confs: conf is binded to sme_conf[0]
	quad_su3 **sme_conf;
        stout_smear_conf_stack_allocate(&sme_conf,conf,pars->stout_pars.nlevls);
        
        //smear iteratively retaining all the stack
        if(phase_pres) addrem_stagphases_to_lx_conf(sme_conf[0]); //remove the staggered phases
        stout_smear_whole_stack(sme_conf,conf,&(pars->stout_pars));
        
        //compute the force in terms of the most smeared conf
        if(phase_pres) addrem_stagphases_to_lx_conf(sme_conf[pars->stout_pars.nlevls]); //add to most smeared conf
	compute_topological_force_lx_conf_internal(F,sme_conf[pars->stout_pars.nlevls],pars,phase_pres);
	
        //remap the force backward
        stouted_force_remap(F,sme_conf,&(pars->stout_pars));
        if(phase_pres) addrem_stagphases_to_lx_conf(sme_conf[0]); //add back again to the original conf
	
	//now free the stack of confs
        stout_smear_conf_stack_free(&sme_conf,pars->stout_pars.nlevls);
      }
    
    //take TA
    gluonic_force_finish_computation(F,conf,phase_pres);
  }
  THREADABLE_FUNCTION_END
}
