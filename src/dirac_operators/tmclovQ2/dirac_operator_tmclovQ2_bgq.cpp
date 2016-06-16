#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ_bgq.hpp"

namespace nissa
{
  void apply_tmclovQ2_bgq(vir_spincolor *out,vir_oct_su3 *conf,double kappa,vir_clover_term_t *Cl,double mu,vir_spincolor *in)
  {
    //application is bufferized
    apply_tmclovQ_bgq(out,conf,kappa,Cl,+mu,in);
    apply_tmclovQ_bgq(out,conf,kappa,Cl,-mu,out);
  }
  
  void apply_tmclovQ2_m2_bgq(vir_spincolor *out,vir_oct_su3 *conf,double kappa,vir_clover_term_t *Cl,double m2,vir_spincolor *in)
  {apply_tmclovQ2_bgq(out,conf,kappa,Cl,sqrt(m2),in);}
  
}
