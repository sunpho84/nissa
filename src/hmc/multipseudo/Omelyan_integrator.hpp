#ifndef _ROOTST_EOIMPR_OMELYAN_INTEGRATOR_HPP
#define _ROOTST_EOIMPR_OMELYAN_INTEGRATOR_HPP

namespace nissa
{
  void evolve_momenta_with_full_rootst_eoimpr_force(quad_su3 **H,quad_su3 **conf,color ***pf,theory_pars_t *theory_pars,double residue,double dt,int *npfs,hmc_force_piece force_piece=BOTH_FORCE_PIECES);
  void omelyan_rootst_eoimpr_evolver(quad_su3 **H,quad_su3 **conf,color ***pf,theory_pars_t *theory_pars,hmc_evol_pars_t *simul);
}

#endif
