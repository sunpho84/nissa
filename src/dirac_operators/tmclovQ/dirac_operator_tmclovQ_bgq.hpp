#ifndef _DIRAC_OPERATOR_TMCLOVQ_BGQ_H
#define _DIRAC_OPERATOR_TMCLOVQ_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_tmclovQ_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double mu,bi_spincolor *in);
}

#endif
