#ifndef _PURE_GAUGE_OMELYAN_INTEGRATOR_HPP
#define _PURE_GAUGE_OMELYAN_INTEGRATOR_HPP

#include "hmc/theory_pars.hpp"

namespace nissa
{
  //parameters for pure gauge theory
  struct pure_gauge_evol_pars_t
  {
    //wheter to use or not hmc
    int use_hmc;
    //basic hmc pars
    double traj_length;
    int skip_mtest_ntraj;
    int nmd_steps;
    //acceleration parameters
    int use_Facc;
    double kappa;
    double residue;
    int naux_fields;
    //number of hb sweeps and hits per link
    int nhb_sweeps;
    int nhb_hits;
    //the same for overrelax
    int nov_sweeps;
    int nov_hits;
    
    pure_gauge_evol_pars_t() : use_hmc(0),traj_length(1.0),skip_mtest_ntraj(30),nmd_steps(13),use_Facc(0),kappa(0.0),residue(1e-12),naux_fields(NDIM),nhb_sweeps(1),nhb_hits(1),nov_sweeps(3),nov_hits(3) {}
  };
  
  void evolve_momenta_with_pure_gauge_force(quad_su3 *H,quad_su3 *conf,theory_pars_t *theory_pars,double dt,quad_su3 *ext_F=NULL);
  void Omelyan_pure_gauge_evolver(quad_su3 *H,quad_su3 *conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *simul);
  void Omelyan_pure_gauge_FACC_evolver(quad_su3 *H,quad_su3 *conf,su3 **pi,su3 **phi,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *simul);
}

#endif
