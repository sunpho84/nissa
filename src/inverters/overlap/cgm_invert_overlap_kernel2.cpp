#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/overlap/dirac_operator_overlap_kernel2.hpp"
#include "linalgs/linalgs.hpp"
#include "inverters/overlap/cgm_invert_overlap_kernel2.hpp"

#define BASETYPE spincolor
#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_overlap_kernel2
#define CGM_OPERATOR_PARAMETERS conf,M,t,

#define CGM_INVERT inv_overlap_kernel2_cgm
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_overlap_kernel_2_cgm_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_overlap_kernel2_cgm

#define CGM_START_COMMUNICATING_BORDERS(A) start_communicating_lx_spincolor_borders(A)
#define CGM_FINISH_COMMUNICATING_BORDERS(A) finish_communicating_lx_spincolor_borders(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION() BASETYPE *t=nissa_malloc("DD_temp",BULK_VOL+BORD_VOL,BASETYPE);
#define CGM_ADDITIONAL_VECTORS_FREE() nissa_free(t);

//additional parameters
#define CGM_NARG 2
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 M
#define CGM_ADDITIONAL_PARAMETERS_CALL conf,M,

#include "inverters/templates/cgm_invert_template_threaded.cpp"
