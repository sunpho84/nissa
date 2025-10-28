#ifndef _SYMANZIK_FORCE_HPP
#define _SYMANZIK_FORCE_HPP

#include <base/field.hpp>

namespace nissa
{
  void Symanzik_force_lx_conf(LxField<quad_su3>& out,
			      const LxField<quad_su3>& conf,
			      const double& beta,
			      const double& C1);
}

#endif
