#ifndef _TMCLOV_FRONTENDS_HPP
#define _TMCLOV_FRONTENDS_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void compute_su3spinspin_tmclov_propagators_multi_mass(su3spinspin ***prop,quad_su3 *conf,double kappa,clover_term_t *Cl,double *mass,int nmass,int niter_max,double *req_res,su3spinspin *source);
}

#endif
