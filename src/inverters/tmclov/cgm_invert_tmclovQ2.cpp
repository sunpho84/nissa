#include <math.h>

#include "../../base/global_variables.h"
#include "../../communicate/communicate.h"
#include "../../base/vectors.h"
#include "../../base/debug.h"
#include "../../dirac_operators/dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define basetype spincolor
#define ndoubles_per_site 24
#define bulk_vol loc_vol
#define bord_vol bord_vol

#define apply_operator apply_tmclovQ2_m2
#define cgm_operator_parameters conf,kappa,csw,Pmunu,t

#define summ_src_and_all_inv_cgm summ_src_and_all_inv_tmclovQ2_m2_cgm
#define cgm_invert inv_tmclovQ2_m2_cgm
#define cgm_invert_run_hm_up_to_comm_prec inv_tmclovQ2_m2_cgm_run_hm_up_to_comm_prec
#define cgm_npossible_requests 16

#define cgm_start_communicating_borders buffered_start_communicating_lx_spincolor_borders
#define cgm_finish_communicating_borders buffered_finish_communicating_lx_spincolor_borders

#define cgm_additional_vectors_allocation()				\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cgm_additional_vectors_free()		\
  nissa_free(t);
#define cgm_additional_parameters_proto quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu
#define cgm_additional_parameters_call conf,kappa,csw,Pmunu

#include "../templates/cgm_invert_template.cpp"

void inv_tmclovQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
{
  double m2[nmass];
  for(int imass=0;imass<nmass;imass++) m2[imass]=m[imass]*m[imass];
  inv_tmclovQ2_m2_cgm(sol,conf,kappa,csw,Pmunu,m2,nmass,niter_max,req_res,source);
}
