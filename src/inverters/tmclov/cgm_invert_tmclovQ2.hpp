#ifndef _CGM_INVERT_CLOVTMQ2_HPP
#define _CGM_INVERT_CLOVTMQ2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
  void inv_tmclovQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
  void inv_tmclovDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
}

#endif
