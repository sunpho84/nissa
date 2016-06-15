#ifndef _DIRAC_OPERATOR_WCLOVQ2_HPP
#define _DIRAC_OPERATOR_WCLOVQ2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_WclovQ2(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor *ext_temp,spincolor *in);
}

#endif
