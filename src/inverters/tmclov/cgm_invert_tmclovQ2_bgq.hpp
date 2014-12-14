#ifndef _CGM_INVERT_TMCLOVQ2_BGQ_H
#define _CGM_INVERT_TMCLOVQ2_BGQ_H

namespace nissa
{
  void inv_tmclovQ2_m2_cgm_bgq(bi_spincolor **sol,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double *m2,int nmass,int niter_max,double *req_res,bi_spincolor *source);
}

#endif
