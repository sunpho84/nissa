#pragma once

#define basetype color
#define ndoubles_per_site 6
#define bulk_vol loc_volh
#define bord_vol bord_volh

#define apply_offdiagonal_operator apply_stD2ee_zero_mass
#define apply_full_operator apply_stD2ee

#define cgmm2s_invert inv_stD2ee_cgmm2s
#define cgmm2s_invert_run_hm_up_to_mach_prec inv_stD2ee_cgmm2s_run_hm_up_to_mach_prec
#define summ_src_and_all_inv_cgmm2s summ_src_and_all_inv_stD2ee_cgmm2s
#define cgmm2s_npossible_requests 16

#define cgmm2s_start_communicating_borders start_communicating_ev_color_borders
#define cgmm2s_finish_communicating_borders finish_communicating_ev_color_borders

#define cgmm2s_additional_vectors_allocation()\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cgmm2s_additional_vectors_free()	\
  nissa_free(t);
#define cgmm2s_additional_parameters_proto quad_su3 **conf,
#define cgmm2s_additional_parameters_call conf,
#define cgmm2s_additional_offdiagonal_parameters conf,t,
#define cgmm2s_additional_full_parameters conf,t,

#include "cgmm2s_invert_template.c"
