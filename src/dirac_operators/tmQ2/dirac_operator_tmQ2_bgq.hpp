#ifndef _DIRAC_OPERATOR_TMQ2_BGQ_HPP
#define _DIRAC_OPERATOR_TMQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmQ2_bgq(vir_spincolor *out,vir_oct_su3 *conf,double kappa,double mu,vir_spincolor *in);
}

#endif
