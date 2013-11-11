#ifndef _DIRAC_OPERATOR_TMQ2_128_BGQ_H
#define _DIRAC_OPERATOR_TMQ2_128_BGQ_H

namespace nissa
{
  void apply_tmQ2_RL_128_bgq(bi_spincolor_128 *out,bi_oct_su3 *conf,double kappa,int RL,double mu,bi_spincolor_128 *in);
}

#endif
