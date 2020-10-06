#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#define BASETYPE vir_spincolor
#define NDOUBLES_PER_SITE 48
#define BULK_VOL loc_volh/2
#define BORD_VOL 0

#define APPLY_OPERATOR tmDkern_eoprec_square_eos_bgq
#define CG_OPERATOR_PARAMETERS temp,conf,kappa,mass,

#define CG_INVERT inv_tmDkern_eoprec_square_eos_cg_64_bgq
#define CG_NPOSSIBLE_REQUESTS 0

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION() \
  vir_spincolor *temp=nissa_malloc("temp",loc_volh/2,vir_spincolor);
#define CG_ADDITIONAL_VECTORS_FREE() \
  nissa_free(temp);

//additional parameters
#define CG_NARG 3
#define AT1 vir_oct_su3**
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 mass

#include "inverters/templates/cg_invert_template_threaded.cpp"
