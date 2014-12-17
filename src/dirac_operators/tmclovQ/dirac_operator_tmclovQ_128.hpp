#ifndef _DIRAC_OPERATOR_TMCLOVQ_128_H
#define _DIRAC_OPERATOR_TMCLOVQ_128_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_tmclovQ_128(spincolor_128 *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,spincolor_128 *in);
  void apply_tmclovQ_128(spincolor_128 *out,quad_su3 *conf,double kappa,opt_as2t_su3 *Cl,double mu,spincolor_128 *in);
}

#endif
