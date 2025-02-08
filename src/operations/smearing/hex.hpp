#ifndef _HEX_HPP
#define _HEX_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void hex_smear_conf(LxField<quad_su3>& sm_conf,
		      const LxField<quad_su3>& conf,
		      const double& alpha1,
		      const double& alpha2,
		      const double& alpha3);
}

#endif
