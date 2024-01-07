#ifndef _MULTIGRIDPARAMS_HPP
#define _MULTIGRIDPARAMS_HPP

#include <cmath>

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
  template <typename T,
	    size_t N>
  struct my_array
  {
    T data[N];
    
    CUDA_HOST_AND_DEVICE inline T& operator[](const size_t i)
    {
      return data[i];
    }
    
    CUDA_HOST_AND_DEVICE const inline T& operator[](const size_t i) const 
    {
      return data[i];
    }
  };
  
  typedef my_array<bool,NDIM> which_dir_t;
  typedef my_array<int,NDIM> coords_t;
  typedef my_array<double,NDIM> momentum_t;
  void read_DDalphaAMG_pars();
  
  namespace multiGrid
  {
    EXTERN_MULTIGRID int nlevels INIT_TO(1);
    EXTERN_MULTIGRID int nsetups[MAX_MG_LEVELS];
    EXTERN_MULTIGRID int smoother_iterations;
    EXTERN_MULTIGRID double mu_factor[MAX_MG_LEVELS];
    EXTERN_MULTIGRID int nu_pre[MAX_MG_LEVELS];
    EXTERN_MULTIGRID int nu_post[MAX_MG_LEVELS];
    EXTERN_MULTIGRID double max_mass INIT_TO(1e300);
    EXTERN_MULTIGRID double max_mass_for_deflation INIT_TO(1e300);
    EXTERN_MULTIGRID bool block_size_set INIT_TO(false);
    EXTERN_MULTIGRID nissa::coords_t block_size[MAX_MG_LEVELS];
    
    EXTERN_MULTIGRID bool setup_valid INIT_TO(false);
    EXTERN_MULTIGRID int use_multiGrid INIT_TO(1);
    EXTERN_MULTIGRID int use_deflated_solver INIT_TO(0);
    EXTERN_MULTIGRID int nEigenvectors INIT_TO(100);
    EXTERN_MULTIGRID double eig_min INIT_TO(6e-2);
    EXTERN_MULTIGRID double eig_max INIT_TO(8);
    
    /// If DDalphaamg is available, check if requested and if the mass is below the maximal
    inline bool checkIfMultiGridAvailableAndRequired(const double& mass)
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
