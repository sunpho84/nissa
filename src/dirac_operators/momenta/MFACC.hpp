#ifndef _OPERATOR_MFACC_H
#define _OPERATOR_MFACC_H

namespace nissa
{
  void apply_MFACC(quad_su3 *out,quad_su3 *conf,double kappa,quad_su3 *in);
  void apply_MFACC(su3 *out,quad_su3 *conf,double kappa,su3 *in);
}

#endif
