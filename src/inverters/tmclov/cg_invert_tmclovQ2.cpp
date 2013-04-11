#include <math.h>

#include "../../base/global_variables.h"
#include "../../base/debug.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define BASETYPE spincolor

#define NDOUBLES_PER_SITE 24
#define BULK_VOL loc_vol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_tmclovQ2
#define CG_OPERATOR_PARAMETERS conf,kappa,csw,Pmunu,temp,mu,

#define CG_INVERT inv_tmclovQ2_cg
#define CG_NPOSSIBLE_REQUESTS 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_lx_spincolor_borders
//#define cg_finish_communicating_borders finish_communicating_lx_spincolor_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()				\
  BASETYPE *temp=nissa_malloc("temp",BULK_VOL+BORD_VOL,BASETYPE);

#define CG_ADDITIONAL_VECTORS_FREE()		\
  nissa_free(temp);

//additional parameters
#define CG_NARG 5
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 csw
#define AT4 as2t_su3*
#define A4 Pmunu
#define AT5 double
#define A5 mu

#include "../templates/cg_invert_template_threaded.cpp"
