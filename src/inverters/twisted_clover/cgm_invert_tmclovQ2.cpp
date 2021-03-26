#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <cmath>

#include "base/vectors.hpp"
#include "base/debug.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor

#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVol
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
    
    inv_tmclovQ2_m2_cgm_portable(sol,conf,kappa,Cl,m2,nmass,niter_max,req_res,source);
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
