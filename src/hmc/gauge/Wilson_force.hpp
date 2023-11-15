#ifndef _WILSON_FORCE_HPP
#define _WILSON_FORCE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/old_field.hpp>
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void Wilson_force_eo_conf(eo_ptr<quad_su3> F,eo_ptr<quad_su3> eo_conf,double beta);
  void Wilson_force_lx_conf(LxField<quad_su3>& F,
			    const LxField<quad_su3>& conf,
			    const double beta);
}

#endif
