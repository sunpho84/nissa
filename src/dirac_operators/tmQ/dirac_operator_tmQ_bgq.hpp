#ifndef _DIRAC_OPERATOR_TMQ_BGQ_H
#define _DIRAC_OPERATOR_TMQ_BGQ_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void hopping_matrix_lx_expand_to_Q_bgq(bi_spincolor *out);
  void hopping_matrix_tmQ_diag_term_bgq(bi_spincolor *out,double kappa,double mu,bi_spincolor *in);
  void apply_tmQ_bgq(bi_spincolor *out,bi_oct_su3 *conf,double kappa,double mu,bi_spincolor *in);
}

#endif
