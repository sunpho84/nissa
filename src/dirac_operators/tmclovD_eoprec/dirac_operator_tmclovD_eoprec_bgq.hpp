#ifndef _DIRAC_OPERATOR_TMCLOV_DEOPREC_BGQ_HPP
#define _DIRAC_OPERATOR_TMCLOV_DEOPREC_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void minus_one_quarter_subt_tmclovDee_or_oo_eos(vir_spincolor *out,double kappa,vir_clover_term_t *Cl,bool dag,double mass,vir_spincolor *in);
  void inv_tmclovDee_or_oo_eos(vir_spincolor *out,vir_inv_clover_term_t *invCl,bool dag,vir_spincolor *in);
  //non squared
  void tmclovDkern_eoprec_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,vir_clover_term_t *Cl_odd,vir_inv_clover_term_t *invCl_evn,bool dag,double mu,vir_spincolor *in);
  
  //squared
  inline void tmclovDkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp1,vir_spincolor *temp2,vir_oct_su3 **conf,double kappa,vir_clover_term_t *Cl_odd,vir_inv_clover_term_t *invCl_evn,double mu,vir_spincolor *in)
  {
    tmclovDkern_eoprec_eos_bgq(temp1,temp2, conf,kappa,Cl_odd,invCl_evn,true, mu,in  );
    tmclovDkern_eoprec_eos_bgq(out,  temp2,conf,kappa,Cl_odd,invCl_evn,false,mu,temp1);
  }
}
#endif
