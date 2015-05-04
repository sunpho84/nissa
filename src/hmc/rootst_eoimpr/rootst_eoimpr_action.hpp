#ifndef _ROOTST_EOIMPR_ACTION_HPP
#define _ROOTST_EOIMPR_ACTION_HPP

namespace nissa
{
  void full_rootst_eoimpr_action(double *res,quad_su3 **eo_conf,quad_su3 **sme_conf,quad_su3 **H,color ***pf,theory_pars_t *theory_pars,hmc_evol_pars_t *simul_pars);
}

#endif
