#ifndef _DIRAC_OPERATOR_TMCLOVQ2_BGQ_H
#define _DIRAC_OPERATOR_TMCLOVQ2_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_tmclovQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *C,double mu,bi_spincolor *in);
  void apply_tmclovQ2_m2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *C,double mu,bi_spincolor *in);    
}

#endif
