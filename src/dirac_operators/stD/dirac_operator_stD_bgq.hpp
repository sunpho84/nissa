#ifndef DIRAC_OPERATOR_STD_BGQ_HPP
#define DIRAC_OPERATOR_STD_BGQ_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_stD2ee_m2_bgq(bi_color* out,bi_oct_su3 **conf,bi_color *temp,double mass2,bi_color *in);
  void apply_single_stD2ee_m2_bgq(bi_single_color* out,bi_single_oct_su3 **conf,bi_single_color *temp,float mass2,bi_single_color *in);
}

#endif
