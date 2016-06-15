#ifndef _DIRAC_OPERATOR_TMCLOVQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ(spincolor *out,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,spincolor *in);
}

#endif
