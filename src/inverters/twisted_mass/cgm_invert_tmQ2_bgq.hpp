#ifndef _CGM_INVERT_TMQ2_BGQ_HPP
#define _CGM_INVERT_TMQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmQ2_m2_cgm_bgq(vir_spincolor **sol,vir_oct_su3 *conf,double kappa,double *m2,int nmass,int niter_max,double *req_res,vir_spincolor *source);
}

#endif
