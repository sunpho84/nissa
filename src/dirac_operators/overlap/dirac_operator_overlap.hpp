#ifndef _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP
#define _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP

#include "new_types/su3.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void rat_approx_for_overlap(quad_su3 *conf, rat_approx_t *appr,double mass_overlap,double maxerr);
  void apply_overlap(spincolor* out,quad_su3 *conf, rat_approx_t *appr,double maxerr,double mass_overlap,double mass, spincolor *in);
}

#endif
