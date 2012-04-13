#pragma once

#define basetype spincolor
#define ndoubles_per_site 24
#define bulk_vol loc_vol
#define bord_vol bord_vol

#define apply_offdiagonal_operator apply_tmQ2_RL
#define apply_full_operator apply_tmQ2_RL

#define summ_src_and_all_inv_cgmm2s summ_src_and_all_inv_tmQ2_RL_cgmm2s
#define cgmm2s_invert inv_tmQ2_RL_cgmm2s
#define cgmm2s_invert_run_hm_up_to_mach_prec inv_tmQ2_RL_cgmm2s_run_hm_up_to_mach_prec
#define cgmm2s_npossible_requests 16

#define cgmm2s_start_communicating_borders start_communicating_lx_spincolor_borders
#define cgmm2s_finish_communicating_borders finish_communicating_lx_spincolor_borders

#define cgmm2s_additional_vectors_allocation()\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cgmm2s_additional_vectors_free()        \
  nissa_free(t);
#define cgmm2s_additional_parameters_proto quad_su3 *conf,double kappa,int RL,
#define cgmm2s_additional_parameters_call conf,kappa,RL,
#define cgmm2s_additional_offdiagonal_parameters conf,kappa,t,RL,0,
#define cgmm2s_additional_full_parameters conf,kappa,t,RL,

#include "cgmm2s_invert_template.c"
