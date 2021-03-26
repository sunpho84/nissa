#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/DDalphaAMG_bridge.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#define BASETYPE spincolor
#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVolh
#define BORD_VOL bord_volh

#define APPLY_OPERATOR tmclovDkern_eoprec_square_eos
#define CG_OPERATOR_PARAMETERS temp1,temp2,eo_conf,kappa,Cl_odd,invCl_evn,mass,

#define CG_INVERT inv_tmclovDkern_eoprec_square_eos_cg_64_portable
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
#define CG_NARG 5
#define AT1 eo_ptr<quad_su3>
#define A1 eo_conf
#define AT2 double
#define A2 kappa
#define AT3 clover_term_t*
#define A3 Cl_odd
#define AT4 inv_clover_term_t*
#define A4 invCl_evn
#define AT5 double
#define A5 mass

#include "inverters/templates/cg_invert_template_threaded.cpp"

namespace nissa
{
  //wrapper
  void inv_tmclovDkern_eoprec_square_eos_cg_64(spincolor *sol,spincolor *guess,eo_ptr<quad_su3> eo_conf,double kappa,double cSW,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,int niter,double residue,spincolor *source)
  {
#if defined USE_DDALPHAAMG
    // if(use_DD and fabs(mu)<=DD::max_mass)
    // 	{
    // 	  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol,quad_su3);
    // 	  paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    // 	  spincolor *tmp_in=nissa_malloc("tmp_in",loc_vol,spincolor);
    // 	  spincolor *tmp_out=nissa_malloc("tmp_out",loc_vol,spincolor);
    // 	  vector_reset(tmp_in);
    // 	  double_vector_copy((double*)tmp_in,(double*)source,loc_volh*sizeof(spincolor)/sizeof(double));
    // 	  DD::solve(tmp_out,lx_conf,kappa,cSW,mu,residue,tmp_in,true);
    // 	  nissa_free(lx_conf);
    // 	  inv_tmclovDkern_eoprec_square_eos_cg_64_portable(sol,guess,eo_conf,kappa,Cl_odd,invCl_evn,mu,niter,residue,source);
    // 	  master_printf("%lg %lg\n",tmp_out[0][0][0][0],sol[0][0][0][0]);
    // 	  master_printf("%lg %lg\n",tmp_out[0][0][0][1],sol[0][0][0][1]);
    // 	  nissa_free(tmp_out);
    // 	  nissa_free(tmp_in);
    // 	}
    // else
    master_printf("DDalpha not yet working, probably expecting a different layout\n");
#endif
    inv_tmclovDkern_eoprec_square_eos_cg_64_portable(sol,guess,eo_conf,kappa,Cl_odd,invCl_evn,mu,niter,residue,source);
  }
}
