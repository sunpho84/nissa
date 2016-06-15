#ifndef _DIRAC_OPERATOR_TMCLOVQ2_128_BGQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_128_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ2_128_bgq(bi_spincolor_128 *out,bi_oct_su3 *conf,double kappa,bi_clover_term_t *Cl,double mu,bi_spincolor_128 *in);
}

#endif
