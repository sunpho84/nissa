#ifndef _ROOTST_EOIMPR_QUARK_FORCE_H
#define _ROOTST_EOIMPR_QUARK_FORCE_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void compute_rootst_eoimpr_quark_and_magnetic_force(quad_su3 **F,double *F_B,quad_su3 **conf,color **pf,theory_pars_t *theory_pars,rat_approx_t *appr,double residue);
}

#endif
