#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

//#include "cg_128_invert_tmQ2.h"

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/debug.h"
#include "../../communicate/communicate.h"
#include "../../dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define BASETYPE bi_spincolor
#define NDOUBLES_PER_SITE 48
#define BULK_VOL loc_volh
#define BORD_VOL 0

#define APPLY_OPERATOR apply_tmQ2_bgq
#define CGM_OPERATOR_PARAMETERS conf,kappa,

#define CGM_INVERT inv_tmQ2_cgm_bgq
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_tmQ2_bgq_cgm_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_tmQ2_bgq_cgm

#define CGM_START_COMMUNICATING_BORDERS(A) 
#define CGM_FINISH_COMMUNICATING_BORDERS(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION() 
#define CGM_ADDITIONAL_VECTORS_FREE() 

//additional parameters
#define CGM_NARG 2
#define AT1 bi_oct_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define CGM_ADDITIONAL_PARAMETERS_CALL conf,kappa,

//#define CG_128_INVERT inv_tmQ2_m2_RL_cg_128
//#define CG_128_ADDITIONAL_PARAMETERS_CALL conf,kappa,RL,

#include "../templates/cgm_invert_template_threaded.cpp"
