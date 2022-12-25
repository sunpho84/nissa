#ifndef _DIRAC_OPERATOR_TMQ2_HPP
#define _DIRAC_OPERATOR_TMQ2_HPP

#include <new_types/su3.hpp>

namespace nissa
{
  void apply_tmQ2(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in);
  
  inline void apply_tmQ2_m2(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double m2,spincolor *in)
  {
    apply_tmQ2(out,conf,kappa,temp,sqrt(m2),in);
  }
}

#endif
