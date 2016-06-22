#ifndef _DIRAC_OPERATOR_TMCLOVQ2_128_BGQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_128_BGQ_HPP

#include "new_types/float_128.hpp"

namespace nissa
{
  void apply_tmclovQ2_128_bgq(vir_spincolor_128 *out,vir_oct_su3 *conf,double kappa,vir_clover_term_t *Cl,double mu,vir_spincolor_128 *in);
}

#endif
