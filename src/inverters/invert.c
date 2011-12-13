#pragma once

//follow inclusion below
#include "cgmms_invert_common.c"

#ifdef BGP

#include "cg_invert_bgp.c"
#include "cgmms_invert_bgp.c"

#else

#include "cg_invert_portable.c"
#include "cgmms_invert_portable.c"

#endif

#include "cg_invert_common.c"

void inv_Q2_cgmms(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit)
{inv_Q2_cgmms_RL(sol,source,guess,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,0);}

void inv_Q2_cgmms_left(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit)
{inv_Q2_cgmms_RL(sol,source,guess,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,1);}
