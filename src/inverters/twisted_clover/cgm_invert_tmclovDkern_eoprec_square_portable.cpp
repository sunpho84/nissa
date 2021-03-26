#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "linalgs/linalgs.hpp"

namespace nissa
{
  void apply_tmclovDkern_cgm_kernel(spincolor *out,spincolor *temp1,spincolor *temp2,eo_ptr<quad_su3> conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,double diag_coeff,spincolor *in)
  {
    tmclovDkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,Cl_odd,invCl_evn,mu,in);
    
    if(diag_coeff!=0)
      double_vector_summassign_double_vector_prod_double((double*)out,(double*)in,diag_coeff,sizeof(spincolor)/sizeof(double)*locVolh);
  }
}


#define BASETYPE spincolor
#define NDOUBLES_PER_SITE 24
#define BULK_VOL locVolh
#define BORD_VOL bord_volh

#define APPLY_OPERATOR apply_tmclovDkern_cgm_kernel
#define CGM_OPERATOR_PARAMETERS t1,t2,conf,kappa,Cl_odd,invCl_evn,mu,

#define CGM_INVERT inv_tmclovDkern_eoprec_square_portable
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_tmclovDkern_eoprec_square_portable_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_tmclovDkern_eoprec_square_portable
#define CGM_NPOSSIBLE_REQUESTS 16

#define CGM_START_COMMUNICATING_BORDERS(A) start_communicating_ev_or_od_spincolor_borders(A,EVN)
#define CGM_FINISH_COMMUNICATING_BORDERS(A) finish_communicating_ev_or_od_spincolor_borders(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION() \
  BASETYPE *t1=nissa_malloc("temp1",BULK_VOL+BORD_VOL,BASETYPE);	\
  BASETYPE *t2=nissa_malloc("temp2",BULK_VOL+BORD_VOL,BASETYPE);

#define CGM_ADDITIONAL_VECTORS_FREE()		\
  nissa_free(t1);				\
  nissa_free(t2);

//additional parameters
#define CGM_NARG 5
#define AT1 eo_ptr<quad_su3>
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 clover_term_t *
#define A3 Cl_odd
#define AT4 inv_clover_term_t *
#define A4 invCl_evn
#define AT5 double
#define A5 mu

#define CGM_ADDITIONAL_PARAMETERS_CALL conf,kappa,Cl_odd,invCl_evn,mu,

#include "inverters/templates/cgm_invert_template_threaded.cpp"

