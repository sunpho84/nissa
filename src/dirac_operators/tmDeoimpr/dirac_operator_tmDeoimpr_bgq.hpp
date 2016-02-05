#ifndef _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void tmDkern_eoprec_square_eos_bgq(bi_spincolor *out,bi_spincolor *temp,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in);
}

#endif
