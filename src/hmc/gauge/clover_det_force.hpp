#ifndef _CLOVER_DET_FORCE_HPP
#define _CLOVER_DET_FORCE_HPP

#include "new_types/su3.hpp"
#include "hmc/quark_pars.hpp"

#include <vector>

namespace nissa
{
  void clover_det_force(quad_su3 **F,std::vector<quark_content_t> quark_content,quad_su3 **eo_conf);
}

#endif
