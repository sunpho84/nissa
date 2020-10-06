#ifndef _DIRAC_OPERATOR_TMCLOVQ2_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ2(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor *temp,double mu,spincolor *in);
  void apply_tmclovQ2_m2(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,spincolor *temp,double mu2,spincolor *in);
}

#endif
