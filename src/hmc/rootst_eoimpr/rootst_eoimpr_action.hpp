#ifndef _rootst_eoimpr_action_hpp
#define _rootst_eoimpr_action_hpp

namespace nissa
{
  void full_rootst_eoimpr_action(double *res,quad_su3 **eo_conf,quad_su3 **sme_conf,quad_su3 **H,double *H_B,color **pf,theory_pars_t *theory_pars,rat_approx_t *appr,double residue);
}

#endif
