#ifndef _ROOTTM_CLOV_EOIMPR_QUARK_FORCE_HPP
#define _ROOTTM_CLOV_EOIMPR_QUARK_FORCE_HPP

#include "new_types/su3.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void summ_the_roottm_clov_eoimpr_quark_force(quad_su3 **F,quad_su3 **eo_conf,double kappa,double cSW,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,spincolor *pf,quad_u1 **u1b,rat_approx_t *appr,double residue);
}

#endif
