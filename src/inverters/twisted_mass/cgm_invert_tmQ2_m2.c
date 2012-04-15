#pragma once

#define basetype spincolor
#define ndoubles_per_site 24
#define bulk_vol loc_vol
#define bord_vol bord_vol

#define apply_operator apply_tmQ2_m2_RL
#define cgm_operator_parameters conf,kappa,t,RL,

#define summ_src_and_all_inv_cgm summ_src_and_all_inv_tmQ2_m2_RL_cgm
#define cgm_invert inv_tmQ2_m2_RL_cgm
#define cgm_invert_run_hm_up_to_mach_prec inv_tmQ2_m2_RL_cgm_run_hm_up_to_mach_prec
#define cgm_npossible_requests 16

#define cgm_start_communicating_borders start_communicating_lx_spincolor_borders
#define cgm_finish_communicating_borders finish_communicating_lx_spincolor_borders

#define cgm_additional_vectors_allocation()\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cgm_additional_vectors_free()        \
  nissa_free(t);
#define cgm_additional_parameters_proto quad_su3 *conf,double kappa,int RL,
#define cgm_additional_parameters_call conf,kappa,RL,

#include "cgm_invert_template.c"
