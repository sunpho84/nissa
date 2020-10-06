#ifndef _TM_FRONTENDS_H
#define _TM_FRONTENDS_H

namespace nissa
{
  void compute_su3spinspin_tm_propagators_multi_mass(su3spinspin ***prop,quad_su3 *conf,double kappa,double *mass,int nmass,int niter_max,double *req_res,su3spinspin *source);
}

#endif
