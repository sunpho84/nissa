#ifndef _CGM_INVERT_TMCLOVQ2_BGQ_HPP
#define _CGM_INVERT_TMCLOVQ2_BGQ_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmclovQ2_m2_cgm_bgq(bi_spincolor **sol,bi_oct_su3 *conf,double kappa,bi_clover_term_t *Cl,double *m2,int nmass,int niter_max,double *req_res,bi_spincolor *source);
}

#endif
