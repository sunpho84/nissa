#include <math.h>

#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define BASETYPE color
#define NDOUBLES_PER_SITE 6
#define BULK_VOL loc_volh
#define BORD_VOL bord_volh

#define APPLY_OPERATOR apply_stD2ee_m2
#define CG_OPERATOR_PARAMETERS conf,t,m*m,

#define CG_INVERT inv_stD2ee_m2_cg
#define CG_NPOSSIBLE_REQUESTS 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()				\
  BASETYPE *t=nissa_malloc("DD_temp",BULK_VOL+BORD_VOL,BASETYPE);
#define CG_ADDITIONAL_VECTORS_FREE()		\
  nissa_free(t);

//additional parameters
#define CGM_NARG 2
#define AT1 quad_su3**
#define A1 conf
#define AT2 double
#define A2 m

#include "../templates/cg_invert_template_threaded.cpp"
