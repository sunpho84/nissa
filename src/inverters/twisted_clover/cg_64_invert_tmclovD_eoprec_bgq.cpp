#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_bgq.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#define BASETYPE vir_spincolor
#define NDOUBLES_PER_SITE 48
#define BULK_VOL loc_volh/2
#define BORD_VOL 0

#define APPLY_OPERATOR tmclovDkern_eoprec_square_eos_bgq
#define CG_OPERATOR_PARAMETERS temp1,temp2,conf,kappa,Cl_odd,invCl_evn,mass,

#define CG_INVERT inv_tmclovDkern_eoprec_square_eos_cg_64_bgq
#define CG_NPOSSIBLE_REQUESTS 0

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION() \
  vir_spincolor *temp1=nissa_malloc("temp1",loc_volh/2,vir_spincolor);	\
  vir_spincolor *temp2=nissa_malloc("temp2",loc_volh/2,vir_spincolor);
#define CG_ADDITIONAL_VECTORS_FREE() \
  nissa_free(temp1);		     \
  nissa_free(temp2);

//additional parameters
#define CG_NARG 5
#define AT1 vir_oct_su3**
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 vir_clover_term_t*
#define A3 Cl_odd
#define AT4 vir_inv_clover_term_t*
#define A4 invCl_evn
#define AT5 double
#define A5 mass

#include "inverters/templates/cg_invert_template_threaded.cpp"
