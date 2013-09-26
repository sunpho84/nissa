#ifndef _DIRAC_OPERATOR_TMCLOVQ2_H
#define _DIRAC_OPERATOR_TMCLOVQ2_H

namespace nissa
{
  void apply_tmclovQ2(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *temp,double mu,spincolor *in);
  void apply_tmclovQ2_m2(spincolor *out,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,spincolor *temp,double mu2,spincolor *in);
}

#endif
