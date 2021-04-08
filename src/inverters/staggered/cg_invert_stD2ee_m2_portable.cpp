#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE color
#define NDOUBLES_PER_SITE 6
#define BULK_VOL locVolh.nastyConvert()
#define BORD_VOL bord_volh

#define APPLY_OPERATOR apply_stD2ee_m2
#define CG_OPERATOR_PARAMETERS conf,t,m2,

#define CG_INVERT inv_stD2ee_m2_cg_portable
#define CG_NPOSSIBLE_REQUESTS 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()				\
  BASETYPE *t=nissa_malloc("DD_temp",BULK_VOL+BORD_VOL,BASETYPE);
#define CG_ADDITIONAL_VECTORS_FREE()		\
  nissa_free(t);

//additional parameters
#define CG_NARG 2
#define AT1 eo_ptr<quad_su3>
#define A1 conf
#define AT2 double
#define A2 m2

#include "inverters/templates/cg_invert_template_threaded.cpp"
