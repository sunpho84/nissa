#ifndef _ROOTST_EOIMPR_QUARK_FORCE_HPP
#define _ROOTST_EOIMPR_QUARK_FORCE_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void compute_rootst_eoimpr_quark_force(quad_su3 **F,quad_su3 **conf,color ***pf,theory_pars_t *theory_pars,rat_approx_t *appr,int *npfs,double residue);
}

#endif
