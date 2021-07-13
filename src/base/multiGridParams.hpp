#ifndef _MULTIGRIDPARAMS_HPP
#define _MULTIGRIDPARAMS_HPP

#include "geometry/geometry_lx.hpp"

#ifndef EXTERN_MULTIGRID
 #define EXTERN_MULTIGRID extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

#ifndef MAX_MG_LEVELS
# define MAX_MG_LEVELS 4
#endif

namespace nissa
{
  void read_DDalphaAMG_pars();
  
  namespace multiGrid
  {
    EXTERN_MULTIGRID int nlevels INIT_TO(1);
    EXTERN_MULTIGRID int nsetups[MAX_MG_LEVELS];
    EXTERN_MULTIGRID int smoother_iterations;
    EXTERN_MULTIGRID double mu_factor[MAX_MG_LEVELS];
    EXTERN_MULTIGRID double max_mass INIT_TO(1e300);
    EXTERN_MULTIGRID bool block_size_set INIT_TO(false);
    EXTERN_MULTIGRID nissa::coords_t block_size[MAX_MG_LEVELS];
    
    EXTERN_MULTIGRID bool setup_valid INIT_TO(false);
    EXTERN_MULTIGRID int use_multiGrid INIT_TO(1);
    EXTERN_MULTIGRID int use_deflated_solver INIT_TO(0);
    
    /// If DDalphaamg is available, check if requested and if the mass is below the maximal
    inline bool checkIfDDalphaAvailableAndRequired(const double& mass)
    {
#if defined USE_DDALPHAAMG or USE_QUDA
      if(use_multiGrid and fabs(mass)<=max_mass)
	return true;
      else
#endif
	return false;
    }
  }
}

#undef EXTERN_MULTIGRID
#undef INIT_TO

#endif
