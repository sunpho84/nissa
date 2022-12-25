#ifndef _CGM_INVERT_TMQ2_HPP
#define _CGM_INVERT_TMQ2_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
  void inv_tmQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
}

#endif
