#ifndef _DIRAC_OPERATOR_TMCLOVQ2_BGQ_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_tmclovQ2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_clover_term_t *Cl,double mu,bi_spincolor *in);
  void apply_tmclovQ2_m2_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,bi_clover_term_t *Cl,double mu,bi_spincolor *in);
}

#endif
