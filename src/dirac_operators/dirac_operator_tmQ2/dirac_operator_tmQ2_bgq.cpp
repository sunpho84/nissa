#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../new_types/new_types_definitions.h"

void apply_tmQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,double mu,bi_spincolor *in)
{
  //application is bufferized
  apply_tmQ2_bgq(out,conf,kappa,+mu,in);
  apply_tmQ2_bgq(out,conf,kappa,-mu,out);
}
