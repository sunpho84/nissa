#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#define BASETYPE spincolor
#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_tmQ2
#define CG_OPERATOR_PARAMETERS conf,kappa,t,m,

#define CG_INVERT inv_tmQ2_cg_64
#define CG_NPOSSIBLE_REQUESTS 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()                              \
  BASETYPE *t=nissa_malloc("DD_temp",locVol+bord_vol,BASETYPE);
#define CG_ADDITIONAL_VECTORS_FREE()            \
  nissa_free(t);

//additional parameters
#define CG_NARG 3
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 m

#include "inverters/templates/cg_invert_template_threaded.cpp"
