#ifndef _THEORY_ACTION_HPP
#define _THEORY_ACTION_HPP
namespace nissa
{
  void compute_quark_action(double *glb_action,quad_su3 **eo_conf,int nfl,std::vector<quad_u1**> u1b, pseudofermion_t *pf,std::vector<quark_content_t> quark_content,hmc_evol_pars_t *simul_pars);
  void full_theory_action(double *res,quad_su3 **eo_conf,quad_su3 **sme_conf,quad_su3 **H,pseudofermion_t *pf,theory_pars_t *theory_pars,hmc_evol_pars_t *simul_pars,double external_quark_action=-1);
}

#endif
