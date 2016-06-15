#ifndef _RECONSTRUCT_TMCLOV_DOUBLET_HPP
#define _RECONSTRUCT_TMCLOV_DOUBLET_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void reconstruct_tmclov_doublet(spincolor *outminus,spincolor *outplus,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,spincolor *in);
}

#endif
