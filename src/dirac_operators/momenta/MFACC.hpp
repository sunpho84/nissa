#ifndef _MFACC_HPP
#define _MFACC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_MFACC(quad_su3 *out,quad_su3 *conf,double kappa,double offset,quad_su3 *in);
  void apply_MFACC(su3 *out,quad_su3 *conf,double kappa,double offset,su3 *in);
  template <class T> void apply_MFACC(T *out,quad_su3 *conf,double kappa,T *in)
  {apply_MFACC(out,conf,kappa,0,in);}
}

#endif
