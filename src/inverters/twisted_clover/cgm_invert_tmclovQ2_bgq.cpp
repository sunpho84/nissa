#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/vectors.hpp"
#include "base/debug.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.hpp"
#include "inverters/twisted_clover/cg_128_invert_tmclovQ2_bgq.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE vir_spincolor
#define NDOUBLES_PER_SITE 48
#define BULK_VOL loc_volh
#define BORD_VOL 0

#define SQRT_SHIFT

#define APPLY_OPERATOR apply_tmclovQ2_bgq
#define CGM_OPERATOR_PARAMETERS conf,kappa,Cl,

#define CGM_INVERT inv_tmclovQ2_m2_cgm_bgq
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_tmclovQ2_cgm_bgq_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_tmclovQ2_cgm_bgq

#define CGM_START_COMMUNICATING_BORDERS(A)
#define CGM_FINISH_COMMUNICATING_BORDERS(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION()
#define CGM_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CGM_NARG 3
#define AT1 vir_oct_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 vir_clover_term_t*
#define A3 Cl
#define CGM_ADDITIONAL_PARAMETERS_CALL conf,kappa,Cl,

#define CG_128_INVERT inv_tmclovQ2_m2_cg_128_bgq
#define CG_128_ADDITIONAL_PARAMETERS_CALL conf,kappa,Cl,

#include "inverters/templates/cgm_invert_template_threaded.cpp"
