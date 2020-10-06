#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/tmQ/dirac_operator_tmQ_bgq.hpp"

#include "routines/ios.hpp"

#include <math.h>

namespace nissa
{
  void apply_tmQ2_bgq(vir_spincolor *out,vir_oct_su3 *conf,double kappa,double mu,vir_spincolor *in)
  {
    //application is bufferized
    apply_tmQ_bgq(out,conf,kappa,+mu,in);
    apply_tmQ_bgq(out,conf,kappa,-mu,out);
  }
}
