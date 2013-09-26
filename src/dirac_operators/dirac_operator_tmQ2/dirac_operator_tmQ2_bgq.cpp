#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "dirac_operator_tmQ/dirac_operator_tmQ_bgq.hpp"

#include "routines/ios.hpp"

#include <math.h>

namespace nissa
{
  void apply_tmQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,double mu,bi_spincolor *in)
  {
    //application is bufferized
    apply_tmQ_bgq(out,conf,kappa,+mu,in);
    apply_tmQ_bgq(out,conf,kappa,-mu,out);
  }
}
