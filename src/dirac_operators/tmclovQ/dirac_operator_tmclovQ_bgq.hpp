#ifndef _DIRAC_OPERATOR_TMCLOVQ_BGQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ_bgq(vir_spincolor *out,vir_oct_su3 *conf,double kappa,vir_clover_term_t *Cl,double mu,vir_spincolor *in);
}

#endif
