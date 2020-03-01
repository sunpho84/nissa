#include <math.h>
#include <cmath>
#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"

#define BASETYPE su3
#define NDOUBLES_PER_SITE 18
#define BULK_VOL loc_vol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_MFACC
#define CGM_OPERATOR_PARAMETERS conf,kappa,

#define CGM_INVERT inv_MFACC_cgm
#define CGM_INVERT_RUN_HM_UP_TO_COMM_PREC inv_MFACC_cgm_run_hm_up_to_comm_prec
#define SUMM_SRC_AND_ALL_INV_CGM summ_src_and_all_inv_MFACC_cgm
#define CGM_NPOSSIBLE_REQUESTS 16

#define CGM_START_COMMUNICATING_BORDERS(A)
#define CGM_FINISH_COMMUNICATING_BORDERS(A)

#define CGM_ADDITIONAL_VECTORS_ALLOCATION()
#define CGM_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CGM_NARG 2
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa

#define CGM_ADDITIONAL_PARAMETERS_CALL conf,kappa,

#include "inverters/templates/cgm_invert_template_threaded.cpp"

namespace nissa
{
  THREADABLE_FUNCTION_7ARG(summ_src_and_all_inv_MFACC_cgm, quad_su3*,sol, quad_su3*,conf, double,kappa, rat_approx_t*,appr, int,niter_max, double,req_res, quad_su3*,source)
  {
    GET_THREAD_ID();
    
    //temporary allocate
    su3 *temp_source=nissa_malloc("temp_source",loc_vol,su3);
    su3 *temp_sol=nissa_malloc("temp_sol",loc_vol,su3);
    
    //do the four links
    for(int mu=0;mu<NDIM;mu++)
      {
	verbosity_lv2_master_printf("Solving Fourier acceleration kernel for internal link %d/4\n",mu);
	
	//copy out
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_copy(temp_source[ivol],source[ivol][mu]);
	set_borders_invalid(temp_source);
	
	//invert
	summ_src_and_all_inv_MFACC_cgm(temp_sol,conf,kappa,appr,niter_max,req_res,temp_source);
	
	//copy in
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_copy(sol[ivol][mu],temp_sol[ivol]);
      }
    set_borders_invalid(sol);
    
    //free temporary storage
    nissa_free(temp_source);
    nissa_free(temp_sol);
  }
  THREADABLE_FUNCTION_END
}
