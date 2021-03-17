#include <math.h>
#include <cmath>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE su3
#define NDOUBLES_PER_SITE 18
#define BULK_VOL loc_vol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_MFACC
#define CG_OPERATOR_PARAMETERS conf,kappa,

#define CG_INVERT inv_MFACC_cg
#define CG_NPOSSIBLE_REQUESTS 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define CG_ADDITIONAL_VECTORS_ALLOCATION()
#define CG_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CG_NARG 2
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa

#include "inverters/templates/cg_invert_template_threaded.cpp"

#include "new_types/su3_op.hpp"

namespace nissa
{
  THREADABLE_FUNCTION_7ARG(inv_MFACC_cg, quad_su3*,sol, quad_su3*,guess, quad_su3*,conf, double,kappa, int,niter, double,residue, quad_su3*,source)
  {
    GET_THREAD_ID();
    
    //temporary allocate
    su3 *temp_source=nissa_malloc("temp_source",loc_vol,su3);
    su3 *temp_sol=nissa_malloc("temp_sol",loc_vol+bord_vol,su3);
    
    //do the four links
    for(int mu=0;mu<NDIM;mu++)
      {
	verbosity_lv2_master_printf("Solving Fourier acceleration kernel for internal link %d/4\n",mu);
	
	//copy out
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_copy(temp_source[ivol],source[ivol][mu]);
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(temp_source);
	
	//invert
	inv_MFACC_cg(temp_sol,NULL,conf,kappa,niter,residue,temp_source);
	
	//copy in
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_copy(sol[ivol][mu],temp_sol[ivol]);
	NISSA_PARALLEL_LOOP_END;
      }
    set_borders_invalid(sol);
    
    //free temporary storage
    nissa_free(temp_source);
    nissa_free(temp_sol);
  }
  THREADABLE_FUNCTION_END
}
