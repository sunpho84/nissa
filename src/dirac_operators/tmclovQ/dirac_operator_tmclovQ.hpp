#ifndef _DIRAC_OPERATOR_TMCLOVQ_H
#define _DIRAC_OPERATOR_TMCLOVQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_tmclovQ(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,spincolor *in);
}

#endif
