#ifndef _DIRAC_OPERATOR_TMQ2_BGQ_HPP
#define _DIRAC_OPERATOR_TMQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,double mu,bi_spincolor *in);
}

#endif
