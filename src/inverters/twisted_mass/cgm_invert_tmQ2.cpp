#include <math.h>
#include <cmath>

#include "cg_128_invert_tmQ2.hpp"

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor
#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVol.nastyConvert()
#define BORD_VOL bordVol.nastyConvert()

#define APPLY_OPERATOR apply_tmQ2_m2_RL
#define CGM_OPERATOR_PARAMETERS conf,kappa,t,RL,

#define CGM_INVERT inv_tmQ2_m2_RL_cgm
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_tmQ2_m2_RL_cgm_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_tmQ2_m2_RL_cgm

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
#define AT3 int
#define A3 RL
#define CGM_ADDITIONAL_PARAMETERS_CALL conf,kappa,RL,

#define CG_128_INVERT inv_tmQ2_m2_RL_cg_128
#define CG_128_ADDITIONAL_PARAMETERS_CALL conf,kappa,RL,

#include "inverters/templates/cgm_invert_template_threaded.cpp"

namespace nissa
{
  void inv_tmQ2_RL_cgm(spincolor **sol,quad_su3 *conf,double kappa,int RL,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {
    double m2[nmass];
    for(int imass=0;imass<nmass;imass++) m2[imass]=m[imass]*m[imass];
    
    inv_tmQ2_m2_RL_cgm(sol,conf,kappa,RL,m2,nmass,niter_max,req_res,source);
  }
  void inv_tmQ2_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {inv_tmQ2_RL_cgm(sol,conf,kappa,0,m,nmass,niter_max,req_res,source);}
  void inv_tmQ2_left_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {inv_tmQ2_RL_cgm(sol,conf,kappa,1,m,nmass,niter_max,req_res,source);}
  
  //put the g5
  void inv_tmDQ_RL_cgm(spincolor** sol,quad_su3* conf,double kappa,int RL,double* m,int nmass,int niter_max,double* req_res,spincolor* source)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int id1=2;id1<4;id1++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ri=0;ri<2;ri++)
	    source[ivol.nastyConvert()][id1][ic1][ri]*=-1;
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
    inv_tmQ2_RL_cgm(sol,conf,kappa,RL,m,nmass,niter_max,req_res,source);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int id1=2;id1<4;id1++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ri=0;ri<2;ri++)
	    source[ivol.nastyConvert()][id1][ic1][ri]*=-1;
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(source);
  }
  
  //wrap RL
  void inv_tmDQ_cgm(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {inv_tmDQ_RL_cgm(sol,conf,kappa,0,m,nmass,niter_max,req_res,source);}
  void inv_tmDQ_cgm_left(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter_max,double *req_res,spincolor *source)
  {inv_tmDQ_RL_cgm(sol,conf,kappa,1,m,nmass,niter_max,req_res,source);}
}
