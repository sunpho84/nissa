#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ_bgq.hpp"

namespace nissa
{
  void apply_tmclovQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double mu,bi_spincolor *in)
  {
    //application is bufferized
    apply_tmclovQ_bgq(out,conf,kappa,Cl,+mu,in);
    apply_tmclovQ_bgq(out,conf,kappa,Cl,-mu,out);
  }

  void apply_tmclovQ2_m2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double m2,bi_spincolor *in)
  {apply_tmclovQ2_bgq(out,conf,kappa,Cl,sqrt(m2),in);}
  
}
