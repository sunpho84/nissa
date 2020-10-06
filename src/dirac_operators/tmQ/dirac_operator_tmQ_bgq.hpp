#ifndef _DIRAC_OPERATOR_TMQ_BGQ_HPP
#define _DIRAC_OPERATOR_TMQ_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void hopping_matrix_lx_expand_to_Q_bgq(vir_spincolor *out);
  void hopping_matrix_tmQ_diag_term_bgq(vir_spincolor *out,double kappa,double mu,vir_spincolor *in);
  void apply_tmQ_bgq(vir_spincolor *out,vir_oct_su3 *conf,double kappa,double mu,vir_spincolor *in);
}

#endif
