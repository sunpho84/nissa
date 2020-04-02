#ifndef _CLOVER_DET_ACTION_HPP
#define _CLOVER_DET_ACTION_HPP

#include "hmc/quark_pars.hpp"

namespace nissa
{
  void clover_det_action(double *act,std::vector<quark_content_t> quark_content,quad_su3 **eo_conf);
}

#endif
