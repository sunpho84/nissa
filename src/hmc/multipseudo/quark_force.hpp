#ifndef _QUARK_FORCE_HPP
#define _QUARK_FORCE_HPP

namespace nissa
{
  void compute_quark_force_finish_computation(quad_su3 **F,quad_su3 **conf);
  void compute_quark_force_no_stout_remapping(quad_su3 **F,quad_su3 **conf,pseudofermion_t *pf,theory_pars_t *tp,rat_approx_t *appr,int *npfs,double residue);
  void compute_quark_force(quad_su3 **F,quad_su3 **conf,pseudofermion_t *pf,theory_pars_t *physics,rat_approx_t *appr,int *npfs,double residue);

}

#endif
