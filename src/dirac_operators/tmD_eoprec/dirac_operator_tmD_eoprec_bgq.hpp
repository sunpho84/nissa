#ifndef _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void tmD_or_Qkern_eoprec_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in,int D_or_Q);
  
  //non squared
  inline void tmDkern_eoprec_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_eos_bgq(out,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_eos_bgq(out,conf,kappa,mu,in,1);}
  
  //squared
  inline void tmD_or_Qkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in,int D_or_Q)
  {
    tmD_or_Qkern_eoprec_eos_bgq(temp,conf,kappa,-mu, in   ,D_or_Q);
    tmD_or_Qkern_eoprec_eos_bgq(out,  conf,kappa,+mu, temp,D_or_Q);
  }
  inline void tmDkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos_bgq(out,temp,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos_bgq(out,temp,conf,kappa,mu,in,1);}
}

#endif
