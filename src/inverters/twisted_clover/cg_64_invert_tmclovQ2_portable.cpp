#include <math.h>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor

#define NDOUBLES_PER_SITE 24
#define BULK_VOL loc_vol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_tmclovQ2
#define CG_OPERATOR_PARAMETERS conf,kappa,Cl,temp,mu,

#define CG_INVERT inv_tmclovQ2_cg_64_portable

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_lx_spincolor_borders
//#define cg_finish_communicating_borders finish_communicating_lx_spincolor_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()				\
  BASETYPE *temp=nissa_malloc("temp",BULK_VOL+BORD_VOL,BASETYPE);

#define CG_ADDITIONAL_VECTORS_FREE()		\
  nissa_free(temp);

//additional parameters
#define CG_NARG 4
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 clover_term_t*
#define A3 Cl
#define AT4 double
#define A4 mu

#include "inverters/templates/cg_invert_template_threaded.cpp"
