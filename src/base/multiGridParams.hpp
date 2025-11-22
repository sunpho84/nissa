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
    EXTERN_MULTIGRID int gcrNkrylov INIT_TO(24);
    
    EXTERN_MULTIGRID double max_mass_for_deflation INIT_TO(1e300);
    EXTERN_MULTIGRID bool block_size_set INIT_TO(false);
    EXTERN_MULTIGRID nissa::Coords block_size[MAX_MG_LEVELS];
    EXTERN_MULTIGRID double max_mass INIT_TO(1e300);
    
    struct SetupPars
    {
      double mu_factor[MAX_MG_LEVELS]{};
      
      double coarse_solver_tol[MAX_MG_LEVELS]{};
      
      int coarse_solver_maxiter[MAX_MG_LEVELS]{};
      
      double smoother_tol[MAX_MG_LEVELS]{};
      
      int nu_pre[MAX_MG_LEVELS]{};
      
      int nu_post[MAX_MG_LEVELS]{};
      
      double omega[MAX_MG_LEVELS]{};
      
      double reliable_delta{0.01};
      
      double reliable_delta_refinement{0.0001};
      
      SetupPars()
      {
	for(int ilev=0;ilev<nlevels;ilev++)
	  nsetups[ilev]=4;
	
	for(int ilev=0;ilev<nlevels;ilev++)
	  mu_factor[ilev]=1;
	
	for(int ilev=0;ilev<nlevels;ilev++)
	  nu_pre[ilev]=0;
	
	for(int ilev=0;ilev<nlevels;ilev++)
	  nu_post[ilev]=7;
    
	constexpr double def_coarse_solver_tol[MAX_MG_LEVELS]={0.15,0.22,0.46};
	
	for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
	  coarse_solver_tol[ilev]=
	  def_coarse_solver_tol[ilev];
	
	constexpr int def_coarse_solver_maxiter[MAX_MG_LEVELS]={100,100,100};
	for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
	  coarse_solver_maxiter[ilev]=
	  def_coarse_solver_maxiter[ilev];
	
	constexpr double def_smoother_tol[MAX_MG_LEVELS]={0.1,0.1,0.15};
	for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
	  smoother_tol[ilev]=
	  def_smoother_tol[ilev];
	
	constexpr double def_omega[MAX_MG_LEVELS]={0.85,0.85,0.85};
	for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
	  omega[ilev]=
	    def_omega[ilev];
      }
    };
    
    EXTERN_MULTIGRID SetupPars pars;
    
    //if no deflation is used, the suffix "no_defl" is not used
    EXTERN_MULTIGRID SetupPars pars_no_defl;
    
#define PROVIDE_ALIAS(X)					\
    inline auto& X=multiGrid::pars.X;				\
    inline auto& X ## _no_deflation=multiGrid::pars_no_defl.X
    
    namespace internal
    {
      PROVIDE_ALIAS(mu_factor);
      PROVIDE_ALIAS(coarse_solver_tol);
      PROVIDE_ALIAS(coarse_solver_maxiter);
      PROVIDE_ALIAS(smoother_tol);
      PROVIDE_ALIAS(nu_pre);
      PROVIDE_ALIAS(nu_post);
      PROVIDE_ALIAS(omega);
      PROVIDE_ALIAS(reliable_delta);
      PROVIDE_ALIAS(reliable_delta_refinement);
    }
    
#undef PROVIDE_ALIAS
    
    EXTERN_MULTIGRID bool setup_valid INIT_TO(false);
    EXTERN_MULTIGRID double setup_refresh_tol INIT_TO(1e3);
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
