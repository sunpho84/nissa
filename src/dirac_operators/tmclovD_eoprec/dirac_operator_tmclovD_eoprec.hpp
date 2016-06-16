#ifndef _DIRAC_OPERATOR_TMCLOVD_EOPREC_HPP
#define _DIRAC_OPERATOR_TMCLOVD_EOPREC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  void tmclovDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  void tmclovDkern_eoprec_eos(spincolor *out,spincolor *temp,quad_su3** conf,double kappa,double mu,spincolor *in);
  void tmclovDkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,double mu,spincolor *in);
  void tmclovn2Deo_eos(spincolor *out,quad_su3 **conf,spincolor *in);
  void tmclovn2Deo_or_tmclovn2Doe_eos(spincolor *out,quad_su3 **conf,int eooe,spincolor *in);
  void tmclovn2Doe_eos(spincolor *out,quad_su3 **conf,spincolor *in);
}

#endif
