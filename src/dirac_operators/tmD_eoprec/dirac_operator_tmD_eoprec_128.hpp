#ifndef _DIRAC_OPERATOR_TMDEOIMPR_128_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_128_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDee_or_oo_eos_128(spincolor_128 *out,double kappa,double mu,spincolor_128 *in);
  void tmDee_or_oo_eos_128(spincolor_128 *out,double kappa,double mu,spincolor_128 *in);
  void tmDkern_eoprec_eos_128(spincolor_128 *out,spincolor *temp,quad_su3** conf,double kappa,double mu,spincolor_128 *in);
  void tmDkern_eoprec_square_eos_128(spincolor_128 *out,spincolor_128 *temp1,spincolor_128 *temp2,quad_su3 **conf,double kappa,double mu,spincolor_128 *in);
}

#endif
