#ifndef _ROOTST_EOIMPR_QUARK_FORCE_HPP
#define _ROOTST_EOIMPR_QUARK_FORCE_HPP

#include "base/field.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void summ_the_rootst_eoimpr_quark_force(EoField<quad_su3>& F,
					  EoField<quad_su3>& eo_conf,
					  const EvnField<color>& pf,
					  const EoField<quad_u1>& u1b,
					  const rat_approx_t& appr,
					  const double& residue);
}

#endif
