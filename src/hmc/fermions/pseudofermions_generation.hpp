#ifndef _PSEUDOFERMIONS_GENERATION_HPP
#define _PSEUDOFERMIONS_GENERATION_HPP

#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"

namespace nissa
{
  void generate_pseudo_fermion(double *action,pseudofermion_t *pf,eo_ptr<quad_su3> conf,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,eo_ptr<quad_u1> u1b,rat_approx_t *rat,double residue,quark_content_t q);
  double generate_pseudofermions(std::vector<std::vector<pseudofermion_t> > &pf,eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,hmc_evol_pars_t &simul_pars,std::vector<rat_approx_t> &rat_appr);
}

#endif
