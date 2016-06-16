#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "communicate/communicate.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/stD/dirac_operator_stD_bgq.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE vir_color
#define NDOUBLES_PER_SITE 12
#define BULK_VOL loc_volh/2
#define BORD_VOL 0

#define APPLY_OPERATOR apply_stD2ee_m2_bgq
#define CGM_OPERATOR_PARAMETERS conf,

#define CGM_INVERT inv_stD2ee_m2_cgm_bgq
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_stD2ee_m2_cgm_bgq_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_stD2ee_m2_cgm_bgq
#define CGM_NPOSSIBLE_REQUESTS 0

#define CGM_START_COMMUNICATING_BORDERS(A)
#define CGM_FINISH_COMMUNICATING_BORDERS(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION()
#define CGM_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CGM_NARG 1
#define AT1 vir_oct_su3**
#define A1 conf

#define CGM_ADDITIONAL_PARAMETERS_CALL conf,

#include "inverters/templates/cgm_invert_template_threaded.cpp"
