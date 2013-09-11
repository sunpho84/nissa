#include <math.h>

#include "communicate/communicate.h"
#include "base/debug.h"
#include "base/global_variables.h"
#include "base/vectors.h"
#include "dirac_operators/dirac_operator_stD/dirac_operator_stD_bgq.h"
#include "linalgs/linalgs.h"
#include "new_types/new_types_definitions.h"

#define BASETYPE bi_color
#define NDOUBLES_PER_SITE 12
#define BULK_VOL loc_volh/2
#define BORD_VOL 0

#define APPLY_OPERATOR apply_stD2ee_m2_bgq
#define CGM_OPERATOR_PARAMETERS conf,t,

#define CGM_INVERT inv_stD2ee_m2_cgm_bgq
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_stD2ee_m2_cgm_bgq_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_stD2ee_m2_cgm_bgq
#define CGM_NPOSSIBLE_REQUESTS 0

#define CGM_START_COMMUNICATING_BORDERS(A)
#define CGM_FINISH_COMMUNICATING_BORDERS(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION() BASETYPE *t=nissa_malloc("DD_temp",BULK_VOL+BORD_VOL,BASETYPE);
#define CGM_ADDITIONAL_VECTORS_FREE() nissa_free(t);

//additional parameters
#define CGM_NARG 1
#define AT1 bi_oct_su3**
#define A1 conf

#define CGM_ADDITIONAL_PARAMETERS_CALL conf,

#include "templates/cgm_invert_template_threaded.cpp"
