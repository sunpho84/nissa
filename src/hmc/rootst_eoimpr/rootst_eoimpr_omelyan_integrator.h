#ifndef _rootst_eoimpr_omeAyan_integrator_h
#define _rootst_eoimpr_omeAyan_integrator_h
void eo_conf_unitarize_explicitly_inverting(quad_su3 **conf);
void evolve_conf_with_momenta(quad_su3 **eo_conf,quad_su3 **H,double dt);
void evolve_momenta_with_full_rootst_eoimpr_force(quad_su3 **H,quad_su3 **conf,color **pf,theory_pars_type *theory_pars,rat_approx_type *appr,double residue,double dt,hmc_force_piece force_piece=BOTH_FORCE_PIECES);
void omelyan_rootst_eoimpr_evolver(quad_su3 **H,quad_su3 **conf,color **pf,theory_pars_type *theory_pars,rat_approx_type *appr,hmc_evol_pars_type *simul,multistep_level multilev=MACRO_STEP);
#endif
