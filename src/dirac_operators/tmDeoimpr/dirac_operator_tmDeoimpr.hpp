#ifndef _DIRAC_OPERATOR_TMDEOIMPR_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  void tmDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  void tmDkern_eoprec_eos(spincolor *out,spincolor *temp,quad_su3** conf,double kappa,double mu,spincolor *in);
  void tmn2Deo_eos(spincolor *out,quad_su3 **conf,spincolor *in);
  void tmn2Deo_or_tmn2Doe_eos(spincolor *out,quad_su3 **conf,int eooe,spincolor *in);
  void tmn2Doe_eos(spincolor *out,quad_su3 **conf,spincolor *in);
  
  void tmD_or_Qkern_eoprec_eos(spincolor *out,spincolor *temp,quad_su3** conf,double kappa,double mu,spincolor *in,int D_or_Q);
  
  //non squared
  inline void tmDkern_eoprec_eos(spincolor *out,spincolor *temp,quad_su3 **conf,double kappa,double mu,spincolor *in)
  {tmD_or_Qkern_eoprec_eos(out,temp,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_eos(spincolor *out,spincolor *temp,quad_su3 **conf,double kappa,double mu,spincolor *in)
  {tmD_or_Qkern_eoprec_eos(out,temp,conf,kappa,mu,in,1);}
  
  //squared
  inline void tmD_or_Qkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,double mu,spincolor *in,int D_or_Q)
  {
    tmD_or_Qkern_eoprec_eos(temp1,temp2,conf,kappa,-mu, in   ,D_or_Q);
    tmD_or_Qkern_eoprec_eos(out,  temp2,conf,kappa,+mu, temp1,D_or_Q);
  }
  inline void tmDkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,double mu,spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,double mu,spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,mu,in,1);}
}

#endif
