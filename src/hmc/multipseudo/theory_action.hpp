#ifndef _THEORY_ACTION_HPP
#define _THEORY_ACTION_HPP

#include "multipseudo_rhmc_step.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void compute_quark_action(double *glb_action,eo_ptr<quad_su3> eo_conf,std::vector<eo_ptr<quad_u1>> u1b, std::vector<std::vector<pseudofermion_t> > *pf,std::vector<quark_content_t> quark_content,hmc_evol_pars_t *simul_pars,std::vector<rat_approx_t> *rat_appr);
  void full_theory_action(double *res,eo_ptr<quad_su3> eo_conf,eo_ptr<quad_su3> sme_conf,eo_ptr<quad_su3> H,std::vector<std::vector<pseudofermion_t> > *pf,theory_pars_t *theory_pars,hmc_evol_pars_t *simul_pars,std::vector<rat_approx_t> *rat_appr,double external_quark_action=-1);
}

#endif
