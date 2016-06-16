#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/stD/dirac_operator_stD_bgq.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE vir_color
#define NDOUBLES_PER_SITE 12
#define BULK_VOL loc_volh/2
#define BORD_VOL 0

#define APPLY_OPERATOR apply_stD2ee_m2_bgq
#define CG_OPERATOR_PARAMETERS conf,m2,

#define CG_INVERT inv_stD2ee_m2_cg_bgq
#define BICGSTAB_INVERT inv_stD2ee_m2_bicgstab_bgq
#define CG_NPOSSIBLE_REQUESTS 0

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()
#define CG_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CG_NARG 2
#define AT1 vir_oct_su3**
#define A1 conf
#define AT2 double
#define A2 m2

#include "inverters/templates/bicgstab_invert_template_threaded.cpp"
#include "inverters/templates/cg_invert_template_threaded.cpp"
