#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#define BASETYPE vir_spincolor
#define NDOUBLES_PER_SITE 48
#define BULK_VOL loc_volh
#define BORD_VOL 0

#define APPLY_OPERATOR apply_tmQ2_bgq
#define CG_OPERATOR_PARAMETERS conf,kappa,m,

#define CG_INVERT inv_tmQ2_RL_cg_64_bgq
#define CG_NPOSSIBLE_REQUESTS 0

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()

#define CG_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CG_NARG 4
#define AT1 vir_oct_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 int
#define A3 RL
#define AT4 double
#define A4 m

#include "inverters/templates/cg_invert_template_threaded.cpp"
