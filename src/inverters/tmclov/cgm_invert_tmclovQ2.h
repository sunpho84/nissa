#ifndef _CGM_INVERT_CLOVTMQ2_H
#define _CGM_INVERT_CLOVTMQ2_H

#include "new_types/new_types_definitions.h"

void inv_tmclovQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
void inv_tmclovQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
void inv_tmclovDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source);
#endif
