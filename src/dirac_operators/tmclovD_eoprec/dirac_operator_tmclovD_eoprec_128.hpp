#ifndef _DIRAC_OPERATOR_TMCLOVD_EOPREC_128_HPP
#define _DIRAC_OPERATOR_TMCLOVD_EOPREC_128_HPP

#include "new_types/float_128.hpp"

namespace nissa
{
  void tmclovDkern_eoprec_square_eos_128(spincolor_128 *out,spincolor_128 *temp1,spincolor_128 *temp2,quad_su3 **conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,spincolor_128 *in);
}

#endif
