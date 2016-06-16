#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/vectors.hpp"
#include "base/debug.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#ifdef BGQ
 #include "geometry/geometry_vir.hpp"
 #include "cgm_invert_tmclovQ2_bgq.hpp"
#endif

#define BASETYPE spincolor

#define NDOUBLES_PER_SITE 24
#define BULK_VOL loc_vol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_tmclovQ2_m2
#define CGM_OPERATOR_PARAMETERS conf,kappa,Cl,t,

#define CGM_INVERT inv_tmclovQ2_m2_cgm_portable
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_tmclovQ2_m2_cgm_run_hm_up_to_comm_prec_portable
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_tmclovQ2_m2_cgm_portable

#define CGM_START_COMMUNICATING_BORDERS(A) start_communicating_lx_spincolor_borders(A)
#define CGM_FINISH_COMMUNICATING_BORDERS(A) finish_communicating_lx_spincolor_borders(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION() BASETYPE *t=nissa_malloc("DD_temp",BULK_VOL+BORD_VOL,BASETYPE);
#define CGM_ADDITIONAL_VECTORS_FREE() nissa_free(t);

//additional parameters
#define CGM_NARG 3
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 clover_term_t*
#define A3 Cl
#define CGM_ADDITIONAL_PARAMETERS_CALL conf,kappa,Cl,

#include "inverters/templates/cgm_invert_template_threaded.cpp"

namespace nissa
{
  void inv_tmclovQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,clover_term_t *Cl,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {
    double m2[nmass];
    for(int imass=0;imass<nmass;imass++) m2[imass]=m[imass]*m[imass];
    
#if defined BGQ
    vir_oct_su3 *vir_conf=nissa_malloc("vir_conf",loc_volh,vir_oct_su3);
    lx_conf_remap_to_virlx(vir_conf,conf);
    vir_clover_term_t *vir_Cl=nissa_malloc("vir_Cl",loc_volh,vir_clover_term_t);
    lx_clover_term_t_remap_to_virlx(vir_Cl,Cl);
    vir_spincolor *vir_source=nissa_malloc("vir_source",loc_volh,vir_spincolor);
    lx_spincolor_remap_to_virlx(vir_source,source);
    vir_spincolor *vir_sol[nmass];
    for(int imass=0;imass<nmass;imass++) vir_sol[imass]=nissa_malloc("vir_sol",loc_volh,vir_spincolor);
    
    inv_tmclovQ2_m2_cgm_bgq(vir_sol,vir_conf,kappa,vir_Cl,m2,nmass,niter_max,req_res,vir_source);
    
    //unmap and free
    for(int imass=0;imass<nmass;imass++)
      {
	virlx_spincolor_remap_to_lx(sol[imass],vir_sol[imass]);
	nissa_free(vir_sol[imass]);
      }
    nissa_free(vir_source);
    nissa_free(vir_conf);
    nissa_free(vir_Cl);
#else
    inv_tmclovQ2_m2_cgm_portable(sol,conf,kappa,Cl,m2,nmass,niter_max,req_res,source);
#endif
  }
  
  void inv_tmclovDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,clover_term_t *Cl,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {
    //put the g5
    NISSA_LOC_VOL_LOOP(ivol) for(int id1=2;id1<NDIRAC;id1++) for(int ic1=0;ic1<NCOL;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
    set_borders_invalid(source);
    inv_tmclovQ2_cgm(sol,conf,kappa,Cl,m,nmass,niter_max,req_res,source);
    NISSA_LOC_VOL_LOOP(ivol) for(int id1=2;id1<NDIRAC;id1++) for(int ic1=0;ic1<NCOL;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
    set_borders_invalid(source);
  }
}
