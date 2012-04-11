#pragma once

#define basetype color
#define ndoubles_per_site 6
#define bulk_vol loc_volh
#define bord_vol bord_volh

#define apply_offdiagonal_operator apply_stD2ee_zero_mass
#define apply_full_operator apply_stD2ee

#define summ_src_and_all_inv_cgmm2s summ_src_and_all_inv_stD2ee_cgmm2s
#define cgmm2s_invert inv_stD2ee_cgmm2s
#define cgmm2s_npossible_requests 16

#define cgmm2s_start_communicating_borders start_communicating_ev_color_borders
#include "cgmms_invert_template.c"
