#ifndef _DIRAC_OPERATOR_TMQ2_128_BGQ_HPP
#define _DIRAC_OPERATOR_TMQ2_128_BGQ_HPP

#include "new_types/float_128.hpp"

namespace nissa
{
  void apply_tmQ2_RL_128_bgq(vir_spincolor_128 *out,vir_oct_su3 *conf,double kappa,int RL,double mu,vir_spincolor_128 *in);
}

#endif
