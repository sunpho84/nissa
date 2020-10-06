#ifndef DIRAC_OPERATOR_STD_BGQ_HPP
#define DIRAC_OPERATOR_STD_BGQ_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_stD2ee_m2_bgq(vir_color* out,vir_oct_su3 **conf,double mass2,vir_color *in);
  void apply_single_stD2ee_m2_bgq(vir_single_color* out,vir_single_oct_su3 **conf,float mass2,vir_single_color *in);
}

#endif
