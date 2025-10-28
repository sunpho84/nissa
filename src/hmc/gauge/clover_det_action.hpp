#ifndef _CLOVER_DET_ACTION_HPP
#define _CLOVER_DET_ACTION_HPP

#include "base/field.hpp"
#include "hmc/quark_pars.hpp"

namespace nissa
{
  double clover_det_action(const std::vector<quark_content_t>& quark_content,
			   const EoField<quad_su3>& conf);
}

#endif
