#ifndef _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP

namespace nissa
{
  void tmD_or_Qkern_eoprec_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in,int D_or_Q);
  
  //non squared
  inline void tmDkern_eoprec_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in)
  {tmD_or_Qkern_eoprec_eos_bgq(out,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in)
  {tmD_or_Qkern_eoprec_eos_bgq(out,conf,kappa,mu,in,1);}
  
  //squared
  inline void tmD_or_Qkern_eoprec_square_eos_bgq(bi_spincolor *out,bi_spincolor *temp,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in,int D_or_Q)
  {
    tmD_or_Qkern_eoprec_eos_bgq(temp,conf,kappa,-mu, in   ,D_or_Q);
    tmD_or_Qkern_eoprec_eos_bgq(out,  conf,kappa,+mu, temp,D_or_Q);
  }
  inline void tmDkern_eoprec_square_eos_bgq(bi_spincolor *out,bi_spincolor *temp,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos_bgq(out,temp,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_square_eos_bgq(bi_spincolor *out,bi_spincolor *temp,bi_oct_su3 **conf,double kappa,double mu,bi_spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos_bgq(out,temp,conf,kappa,mu,in,1);}
}

#endif
