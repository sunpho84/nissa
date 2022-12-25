#ifndef _PURE_GAUGE_OMELYAN_INTEGRATOR_HPP
#define _PURE_GAUGE_OMELYAN_INTEGRATOR_HPP

#include "hmc/theory_pars.hpp"

namespace nissa
{
  //parameters for pure gauge theory
  struct pure_gauge_evol_pars_t
  {
    //whether to use or not hmc
    int use_hmc;
    
    //basic hmc pars
    double traj_length;
    int skip_mtest_ntraj;
    int nmd_steps;
    
    //number of hb sweeps and hits per link
    int nhb_sweeps;
    int nhb_hits;
    
    //the same for overrelax
    int nov_sweeps;
    int nov_hits;
    
    int def_use_hmc(){return 0;}
    int def_traj_length(){return 1.0;}
    int def_skip_mtest_ntraj(){return 30;}
    int def_nmd_steps(){return 11;}
    
    //number of hb sweeps and hits per link
    int def_nhb_sweeps(){return 1;}
    int def_nhb_hits(){return 1;}
    
    //the same for overrelax
    int def_nov_sweeps(){return 3;}
    int def_nov_hits(){return 3;}
    
    pure_gauge_evol_pars_t() : use_hmc(def_use_hmc()),traj_length(def_traj_length()),skip_mtest_ntraj(def_skip_mtest_ntraj()),nmd_steps(def_nmd_steps()),nhb_sweeps(def_nhb_sweeps()),nhb_hits(def_nhb_hits()),nov_sweeps(def_nov_sweeps()),nov_hits(def_nov_hits()) {}
    
    std::string get_str(int full=false)
    {
      std::ostringstream os;
      
      if(full or is_nonstandard())
	{
	  os<<"QuenchedEvolution\n";
	  if(full or use_hmc!=def_use_hmc()) os<<"UseHMC\t=\t"<<use_hmc<<"\n";
	  if(full or traj_length!=def_traj_length()) os<<"TrajLength\t=\t"<<traj_length<<"\n";
	  if(full or skip_mtest_ntraj!=def_skip_mtest_ntraj()) os<<"SkipMetro\t=\t"<<skip_mtest_ntraj<<"\n";
	  if(full or nmd_steps!=def_nmd_steps()) os<<"NSteps\t\t=\t"<<nmd_steps<<"\n";
	  if(full or nhb_sweeps!=def_nhb_sweeps()) os<<"NHBSweeps\t=\t"<<nhb_sweeps<<"\n";
	  if(full or nhb_hits!=def_nhb_hits()) os<<"NHBHits\t\t=\t"<<nhb_hits<<"\n";
	  if(full or nov_sweeps!=def_nov_sweeps()) os<<"NOVSweeps\t=\t"<<nov_sweeps<<"\n";
	  if(full or nov_hits!=def_nov_hits()) os<<"NOVHits\t\t=\t"<<nov_hits<<"\n";
	}
      
      return os.str();
    }
    
    bool is_nonstandard()
    {
      return
	use_hmc!=def_use_hmc() or
	traj_length!=def_traj_length() or
	skip_mtest_ntraj!=def_skip_mtest_ntraj() or
	nmd_steps!=def_nmd_steps() or
	nhb_sweeps!=def_nhb_sweeps() or
	nhb_hits!=def_nhb_hits() or
	nov_sweeps!=def_nov_sweeps() or
	nov_hits!=def_nov_hits();
    }
  };
  
  void evolve_momenta_with_pure_gauge_force(LxField<quad_su3>& H,
					    const LxField<quad_su3>& conf,
					    const theory_pars_t& theory_pars,
					    const double& dt,
					    LxField<quad_su3>& F);
  
  void Omelyan_pure_gauge_evolver(quad_su3 *H,quad_su3 *conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *simul);
}

#endif
