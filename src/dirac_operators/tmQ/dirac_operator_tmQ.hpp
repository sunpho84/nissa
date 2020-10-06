#ifndef _DIRAC_OPERATOR_TMQ_HPP
#define _DIRAC_OPERATOR_TMQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmQ_RL(spincolor *out,quad_su3 *conf,double kappa,double mu,int RL,spincolor *in);
  void apply_tmQ(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in);
  void apply_tmQ_v1(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in);
}

#endif
