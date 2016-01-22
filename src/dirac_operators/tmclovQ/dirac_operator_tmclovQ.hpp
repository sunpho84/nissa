#ifndef _DIRAC_OPERATOR_TMCLOVQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double mu,spincolor *in);
}

#endif
