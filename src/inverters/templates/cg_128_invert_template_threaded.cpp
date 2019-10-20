// Template to invert in quadruple precision using c.g.
// The solution is obtained in a mixed precision approach.

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#if CG_NARG >= 6
 #error not supported
#endif

namespace nissa
{
#if CG_NARG == 0
  THREADABLE_FUNCTION_5ARG(CG_128_INVERT, BASETYPE*,sol, BASETYPE*,guess, int,niter, double,external_solver_residue, BASETYPE*,external_source)
#elif CG_NARG == 1
  THREADABLE_FUNCTION_6ARG(CG_128_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, int,niter, double,external_solver_residue, BASETYPE*,external_source)
#elif CG_NARG == 2
  THREADABLE_FUNCTION_7ARG(CG_128_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, int,niter, double,external_solver_residue, BASETYPE*,external_source)
#elif CG_NARG == 3
  THREADABLE_FUNCTION_8ARG(CG_128_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, AT3,A3, int,niter, double,external_solver_residue, BASETYPE*,external_source)
#elif CG_NARG == 4
  THREADABLE_FUNCTION_9ARG(CG_128_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, AT3,A3, AT4,A4, int,niter, double,external_solver_residue, BASETYPE*,external_source)
#elif CG_NARG == 5
  THREADABLE_FUNCTION_10ARG(CG_128_INVERT, BASETYPE*,sol, BASETYPE*,guess, AT1,A1, AT2,A2, AT3,A3, AT4,A4, AT5,A5, int,niter, double,external_solver_residue, BASETYPE*,external_source)
#endif
  {
    //Allocate the solution in 128 bit and initialize it. If guess passed copy it.
    BASETYPE_128 *sol_128=nissa_malloc("sol_128",BULK_SIZE+BORD_SIZE,BASETYPE_128);
    CG_ADDITIONAL_VECTORS_ALLOCATION();
    if(guess==NULL) vector_reset(sol_128);
    else quadruple_vector_from_double_vector((float_128*)sol_128,(double*)guess,BULK_SIZE*NDOUBLES_PER_SITE);
    
    //compute and print source norm
    double source_norm;
    double_vector_glb_scalar_prod(&source_norm,(double*)external_source,(double*)external_source,BULK_SIZE*NDOUBLES_PER_SITE);
    verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
    
    //internal inverter source
    BASETYPE *internal_source=nissa_malloc("internal_source",BULK_SIZE+BORD_SIZE,BASETYPE);
    
    //residue in 128 bit
    BASETYPE_128 *residue_128=nissa_malloc("residue_128",BULK_SIZE+BORD_SIZE,BASETYPE_128);
    
    //temporary vector for inverter
    BASETYPE_128 *temp_128=nissa_malloc("temp",BULK_SIZE+BORD_SIZE,BASETYPE_128); 
    
    //invert until residue is lower than required
    double previous_residue=1,current_residue;
    int quit=0,ext_iter=0;
    do
      {
	// 1) compute D*sol in quadruple precision
	APPLY_OPERATOR_128(residue_128,CG_OPERATOR_128_PARAMETERS sol_128);
	
	// 2) compute the new internal_source = external_source - OP * sol_128, and residue module
	quadruple_vector_subt_from_double_vector((float_128*)residue_128,(double*)external_source,(float_128*)residue_128,BULK_SIZE*NDOUBLES_PER_SITE);
	
	double_vector_from_quadruple_vector((double*)internal_source,(float_128*)residue_128,BULK_SIZE*NDOUBLES_PER_SITE);
	
	double_conv_quadruple_vector_glb_scalar_prod(&current_residue,(float_128*)residue_128,(float_128*)residue_128,BULK_SIZE*NDOUBLES_PER_SITE);
	current_residue/=source_norm;
	verbosity_lv2_master_printf("\nExternal loop iter %d relative residue: %lg\n\n",ext_iter,current_residue);
	
	// 3) calibrate inner solver stopping condition
	double inner_solver_residue=std::max(1.e-16,external_solver_residue/current_residue);
	
	// 3) if residue not reached, compute the new approximated solution
	if(current_residue>=external_solver_residue)
	  {
	    //compute partial sol
	    CG_128_INNER_SOLVER(sol,NULL,CG_128_INNER_PARAMETERS_CALL 1000000,inner_solver_residue,internal_source);
	    
	    //add the approximated solution to the total one
	    quadruple_vector_summassign_double_vector((float_128*)sol_128,(double*)sol,BULK_SIZE*NDOUBLES_PER_SITE);
	  }
	if(ext_iter!=0 && !(current_residue<previous_residue))
	  {
	    quit=1;
	    master_printf("Previous residue %lg, current residue %lg. Not converging, quitting loop\n",previous_residue,current_residue);
	  }
	previous_residue=current_residue;
	
	ext_iter++;
      }
    while(current_residue>=external_solver_residue && !quit);
    
    verbosity_lv1_master_printf("\nFinal residue: %lg\n",current_residue);
    verbosity_lv2_master_printf("\n");
    
    //copy the solution in 128 bit to the 64 bit
    double_vector_from_quadruple_vector((double*)sol,(float_128*)sol_128,BULK_SIZE*NDOUBLES_PER_SITE);
    
    CG_ADDITIONAL_VECTORS_FREE();
    
    nissa_free(residue_128);
    nissa_free(temp_128);
    nissa_free(sol_128);
    nissa_free(internal_source);
  }
  THREADABLE_FUNCTION_END
}

#undef BASETYPE
#undef BASETYPE_128
#undef NDOUBLES_PER_SITE
#undef BULK_SIZE
#undef BORD_SIZE

#undef APPLY_OPERATOR_128
#undef CG_OPERATOR_128_PARAMETERS
#undef CG_128_INVERT
#undef CG_128_INNER_PARAMETERS_CALL
#undef CG_128_INNER_SOLVER
#undef CG_ADDITIONAL_VECTORS_FREE
#undef CG_ADDITIONAL_VECTORS_ALLOCATE
