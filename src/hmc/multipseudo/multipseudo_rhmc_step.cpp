#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <math.h>
#include <vector>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

#include "hmc/backfield.hpp"
#include "hmc/hmc.hpp"
#include "hmc/fermions/pseudofermions_generation.hpp"
#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/multipseudo/Omelyan_integrator.hpp"
#include "hmc/multipseudo/set_expansions.hpp"
#include "hmc/multipseudo/theory_action.hpp"


namespace nissa
{
  //ROOT_STAG has mass term completely constant
  namespace
  {
    const double shift_poles=1,shift_poles_back=-1;
    void shift_all_ROOT_STAG_poles(theory_pars_t &theory_pars,std::vector<rat_approx_t> &rat_appr,const double sign)
    {
      
      for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
	if(theory_pars.quarks[iflav].discretiz==ferm_discretiz::ROOT_STAG)
	  for(int i=0;i<nappr_per_quark;i++)
	    rat_appr[iflav*nappr_per_quark+i].shift_all_poles(sign*sqr(theory_pars.quarks[iflav].mass));
    }
  }
  
  //perform a full hmc step and return the difference between final and original action
  double multipseudo_rhmc_step(eo_ptr<quad_su3> out_conf,eo_ptr<quad_su3> in_conf,theory_pars_t &theory_pars,hmc_evol_pars_t &simul_pars,std::vector<rat_approx_t> &rat_appr,int itraj)
  {
    //header
    master_printf("Trajectory %d->%d\n",itraj,itraj+1);
    master_printf("-----------------------\n");
    
    //take initial time
    double hmc_time=-take_time();
    
    //allocate the momenta
    eo_ptr<quad_su3> H;
    for(int par=0;par<2;par++) H[par]=nissa_malloc("H",loc_volh,quad_su3);
    
    //copy the old conf into the new
    for(int par=0;par<2;par++)
      {
	vector_copy(out_conf[par],in_conf[par]);
	set_borders_invalid(out_conf[par]);
      }
    
    //allocate pseudo-fermions
    std::vector<std::vector<pseudofermion_t> > pf(theory_pars.nflavs());
    for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
      {
	int npf=simul_pars.npseudo_fs[iflav];
	pf[iflav].resize(npf);
	for(int ipf=0;ipf<npf;ipf++) pf[iflav][ipf].create(theory_pars.quarks[iflav].discretiz);
      }
    
    //if needed smear the configuration for pseudo-fermions, approx generation and action computation
    //otherwise bind out_conf to sme_conf
    eo_ptr<quad_su3> sme_conf;
    for(int eo=0;eo<2;eo++)
      sme_conf[eo]=(theory_pars.stout_pars.nlevels!=0)?
	nissa_malloc("sme_conf",loc_volh+bord_volh+edge_volh,quad_su3):out_conf[eo];
    if(theory_pars.stout_pars.nlevels!=0)
      {
	verbosity_lv2_master_printf("Stouting the links for pseudo-fermions generation and initial action computation\n");
	stout_smear(sme_conf,out_conf,&(theory_pars.stout_pars));
	
	verbosity_lv2_master_printf("Original plaquette: %16.16lg\n",global_plaquette_eo_conf(out_conf));
	verbosity_lv2_master_printf("Stouted plaquette: %16.16lg\n",global_plaquette_eo_conf(sme_conf));
      }
    
    //generate the appropriate expansion of rational approximations
    set_expansions(&rat_appr,sme_conf,&theory_pars,&simul_pars);
    
    //shift all the poles of the mass for staggered operator
    shift_all_ROOT_STAG_poles(theory_pars,rat_appr,shift_poles);
    
    //generate all pseudofermions and momenta
    double pf_action=generate_pseudofermions(pf,sme_conf,theory_pars,simul_pars,rat_appr);
    generate_hmc_momenta(H);
    
    //compute initial action
    double init_action;
    full_theory_action(&init_action,out_conf,sme_conf,H,&pf,&theory_pars,&simul_pars,&rat_appr,pf_action);
    verbosity_lv2_master_printf("Initial action: %lg\n",init_action);
    
    //evolve
    Omelyan_integrator(H,out_conf,&pf,&theory_pars,&simul_pars,&rat_appr);
    
    //if needed, resmear the conf
    if(theory_pars.stout_pars.nlevels!=0)
      {
	verbosity_lv2_master_printf("Stouting the links for final action computation\n");
	stout_smear(sme_conf,out_conf,&(theory_pars.stout_pars));
      }
    
    //compute final action
    double final_action;
    full_theory_action(&final_action,out_conf,sme_conf,H,&pf,&theory_pars,&simul_pars,&rat_appr);
    verbosity_lv2_master_printf("Final action: %lg\n",final_action);
    
    //compute the diff
    double diff_action=final_action-init_action;
    
    //free smeared conf
    for(int par=0;par<2;par++) nissa_free(H[par]);
    if(theory_pars.stout_pars.nlevels!=0) for(int eo=0;eo<2;eo++) nissa_free(sme_conf[eo]);
    
    //shift back all the poles of the mass for staggered operator
    shift_all_ROOT_STAG_poles(theory_pars,rat_appr,shift_poles_back);
    
    //take time
    hmc_time+=take_time();
    verbosity_lv1_master_printf("Total time to perform rhmc step: %lg s\n",hmc_time);
    
    return diff_action;
  }
}
