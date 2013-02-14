#include <math.h>

#include "cg_128_invert_tmQ2.h"

#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../../base/debug.h"
#include "../../dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define basetype spincolor
#define ndoubles_per_site 24
#define bulk_vol loc_vol
#define bord_vol bord_vol

#define apply_operator apply_tmQ2_m2_RL
#define cgm_operator_parameters conf,kappa,t,RL

#define summ_src_and_all_inv_cgm summ_src_and_all_inv_tmQ2_m2_RL_cgm
#define cgm_invert inv_tmQ2_m2_RL_cgm
#define cgm_invert_run_hm_up_to_comm_prec inv_tmQ2_m2_RL_cgm_run_hm_up_to_comm_prec
#define cgm_npossible_requests 16

#define cgm_start_communicating_borders start_communicating_lx_spincolor_borders
#define cgm_finish_communicating_borders finish_communicating_lx_spincolor_borders

#define cgm_additional_vectors_allocation()				\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cgm_additional_vectors_free()		\
  nissa_free(t);
#define cgm_additional_parameters_proto quad_su3 *conf,double kappa,int RL
#define cgm_additional_parameters_call conf,kappa,RL

#define cg_128_invert inv_tmQ2_m2_RL_cg_128
#define cg_128_additional_parameters_call conf,kappa,RL

#include "../templates/cgm_invert_template.cpp"


void inv_tmQ2_RL_cgm(spincolor **sol,quad_su3 *conf,double kappa,int RL,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{
  double m2[nmass];
  for(int imass=0;imass<nmass;imass++) m2[imass]=m[imass]*m[imass];
  inv_tmQ2_m2_RL_cgm(sol,conf,kappa,RL,m2,nmass,niter_max,req_res,source);
}
void inv_tmQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{inv_tmQ2_RL_cgm(sol,conf,kappa,0,m,nmass,niter_max,req_res,source);}
void inv_tmQ2_left_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{inv_tmQ2_RL_cgm(sol,conf,kappa,1,m,nmass,niter_max,req_res,source);}

void inv_tmDQ_RL_cgm(spincolor **sol,quad_su3 *conf,double kappa,int RL,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{
  //put the g5
  nissa_loc_vol_loop(ivol) for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
  set_borders_invalid(source);
  inv_tmQ2_RL_cgm(sol,conf,kappa,RL,m,nmass,niter_max,req_res,source);
  nissa_loc_vol_loop(ivol) for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
  set_borders_invalid(source);
}

void inv_tmDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{inv_tmDQ_RL_cgm(sol,conf,kappa,0,m,nmass,niter_max,req_res,source);}
void inv_tmDQ_cgm_left(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{inv_tmDQ_RL_cgm(sol,conf,kappa,1,m,nmass,niter_max,req_res,source);}

