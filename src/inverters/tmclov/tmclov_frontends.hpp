#ifndef _TMCLOV_FRONTENDS_H
#define _TMCLOV_FRONTENDS_H

namespace nissa
{
  void compute_su3spinspin_tmclov_propagators_multi_mass(su3spinspin ***prop,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *mass,int nmass,int niter_max,double *req_res,su3spinspin *source);
}

#endif
