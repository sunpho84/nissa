#include <math.h>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE single_color
#define NDOUBLES_PER_SITE 6
#define BULK_VOL loc_volh
#define BORD_VOL bord_volh

#define APPLY_OPERATOR apply_stD2ee_m2_32
#define CGM_OPERATOR_PARAMETERS conf,t,

#define CGM_INVERT inv_stD2ee_m2_cgm_32_portable
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_stD2ee_m2_cgm_32_portable_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_stD2ee_m2_cgm_32_portable
#define CGM_NPOSSIBLE_REQUESTS 16

#define CGM_START_COMMUNICATING_BORDERS(A) start_communicating_ev_or_od_single_color_borders(A,EVN)
#define CGM_FINISH_COMMUNICATING_BORDERS(A) finish_communicating_ev_or_od_single_color_borders(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION() BASETYPE *t=nissa_malloc("DD_temp",BULK_VOL+BORD_VOL,BASETYPE);
#define CGM_ADDITIONAL_VECTORS_FREE() nissa_free(t);

//additional parameters
#define CGM_NARG 1
#define AT1 single_quad_su3**
#define A1 conf

#define CGM_ADDITIONAL_PARAMETERS_CALL conf,

#include "inverters/templates/cgm_32_invert_template_threaded.cpp"
