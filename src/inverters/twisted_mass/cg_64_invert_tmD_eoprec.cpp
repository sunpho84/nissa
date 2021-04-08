#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#define BASETYPE spincolor
#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVolh.nastyConvert()
#define BORD_VOL bord_volh

#define APPLY_OPERATOR tmDkern_eoprec_square_eos
#define CG_OPERATOR_PARAMETERS temp1,temp2,conf,kappa,mass,

#define CG_INVERT inv_tmDkern_eoprec_square_eos_cg_64_portable
#define CG_NPOSSIBLE_REQUESTS 0

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()                              \
  BASETYPE *temp1=nissa_malloc("temp1",BULK_VOL+BORD_VOL,BASETYPE); \
  BASETYPE *temp2=nissa_malloc("temp2",BULK_VOL+BORD_VOL,BASETYPE);

#define CG_ADDITIONAL_VECTORS_FREE()            \
  nissa_free(temp1);				\
  nissa_free(temp2);

//additional parameters
#define CG_NARG 3
#define AT1 eo_ptr<quad_su3>
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 mass

#include "inverters/templates/cg_invert_template_threaded.cpp"

namespace nissa
{
  //wrapper for bgq
  void inv_tmDkern_eoprec_square_eos_cg_64(spincolor *sol,spincolor *guess,eo_ptr<quad_su3> eo_conf,double kappa,double mu,int niter,double residue,spincolor *source)
  {
    inv_tmDkern_eoprec_square_eos_cg_64_portable(sol,guess,eo_conf,kappa,mu,niter,residue,source);
  }
}
