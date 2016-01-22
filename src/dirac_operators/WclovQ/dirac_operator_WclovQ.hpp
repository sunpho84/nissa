#ifndef _DIRAC_OPERATOR_WCLOVQ_HPP
#define _DIRAC_OPERATOR_WCLOVQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_WclovQ(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *in);
}

#endif
