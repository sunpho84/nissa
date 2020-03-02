#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <cmath>

#ifdef BGQ
 #include "cg_64_invert_tmD_eoprec_bgq.hpp"
 #include "geometry/geometry_eo.hpp"
 #include "geometry/geometry_vir.hpp"
#endif
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
#define BULK_VOL loc_volh
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
#define AT1 quad_su3**
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 mass

#include "inverters/templates/cg_invert_template_threaded.cpp"

namespace nissa
{
  //wrapper for bgq
  void inv_tmDkern_eoprec_square_eos_cg_64(spincolor *sol,spincolor *guess,quad_su3 **eo_conf,double kappa,double mu,int niter,double residue,spincolor *source)
  {
#ifdef BGQ
    //allocate
    vir_spincolor *vir_source=nissa_malloc("vir_source",loc_volh/2,vir_spincolor);
    vir_oct_su3 *vir_eo_conf[2]={nissa_malloc("vir_conf_evn",loc_volh+bord_volh,vir_oct_su3),
                               nissa_malloc("vir_conf_odd",loc_volh+bord_volh,vir_oct_su3)};
    vir_spincolor *vir_sol=nissa_malloc("vir_sol",loc_volh/2,vir_spincolor);
    vir_spincolor *vir_guess=(guess!=NULL)?nissa_malloc("vir_guess",loc_volh/2,vir_spincolor):NULL;
    
    ////////////////////////
    
    //remap in
    evn_or_odd_spincolor_remap_to_virevn_or_odd(vir_source,source,ODD);
    eo_conf_remap_to_vireo(vir_eo_conf,eo_conf);
    if(guess!=NULL) evn_or_odd_spincolor_remap_to_virevn_or_odd(vir_guess,guess,ODD);
    
    //invert
    inv_tmDkern_eoprec_square_eos_cg_64_bgq(vir_sol,vir_guess,vir_eo_conf,kappa,mu,niter,residue,vir_source);
    
    //remap out
    virevn_or_odd_spincolor_remap_to_evn_or_odd(sol,vir_sol,ODD);
    
    ////////////////////////
    
    //free
    nissa_free(vir_eo_conf[EVN]);
    nissa_free(vir_eo_conf[ODD]);
    nissa_free(vir_source);
    nissa_free(vir_sol);
    if(guess!=NULL) nissa_free(vir_guess);
#else
    inv_tmDkern_eoprec_square_eos_cg_64_portable(sol,guess,eo_conf,kappa,mu,niter,residue,source);
#endif
  } 
}
