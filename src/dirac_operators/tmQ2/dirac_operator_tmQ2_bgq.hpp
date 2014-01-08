#ifndef _DIRAC_OPERATOR_TMQ2_BGQ_H
#define _DIRAC_OPERATOR_TMQ2_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_tmQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,double mu,bi_spincolor *in);
}

#endif
