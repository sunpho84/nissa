#ifndef _DIRAC_OPERATOR_TMD_EOIMPR_HPP
#define _DIRAC_OPERATOR_TMD_EOIMPR_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  void tmDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  void tmDkern_eoprec_eos(spincolor *out,spincolor *temp,quad_su3** conf,double kappa,double mu,spincolor *in);
  void tmDkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,double mu,spincolor *in);
  void tmn2Deo_eos(spincolor *out,quad_su3 **conf,spincolor *in);
  void tmn2Deo_or_tmn2Doe_eos(spincolor *out,quad_su3 **conf,int eooe,spincolor *in);
  void tmDkern_eoprec_eos_put_together_and_include_gamma5(spincolor *out,spincolor *temp);
  void tmn2Doe_eos(spincolor *out,quad_su3 **conf,spincolor *in);
}

#endif
