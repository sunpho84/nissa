#pragma once

#define basetype color
#define ndoubles_per_site 6
#define bulk_vol loc_volh
#define bord_vol bord_volh

#define apply_operator apply_stD2ee_m2
#define cgm_operator_parameters conf,t,

#define cgm_invert inv_stD2ee_m2_cgm
#define cgm_invert_run_hm_up_to_mach_prec inv_stD2ee_m2_cgm_run_hm_up_to_mach_prec
#define summ_src_and_all_inv_cgm summ_src_and_all_inv_stD2ee_m2_cgm
#define cgm_npossible_requests 16

#define cgm_start_communicating_borders start_communicating_ev_color_borders
#define cgm_finish_communicating_borders finish_communicating_ev_color_borders

#define cgm_additional_vectors_allocation()\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cgm_additional_vectors_free()	\
  nissa_free(t);
#define cgm_additional_parameters_proto quad_su3 **conf,
#define cgm_additional_parameters_call conf,

#include "cgm_invert_template.c"
