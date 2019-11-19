#ifndef _ROOTST_EOIMPR_QUARK_FORCE_HPP
#define _ROOTST_EOIMPR_QUARK_FORCE_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void summ_the_rootst_eoimpr_quark_force(eo_ptr<quad_su3> F,double charge,eo_ptr<quad_su3> eo_conf,color *pf,int quantization,eo_ptr<quad_u1> u1b,rat_approx_t *appr,double residue);
}

#endif
