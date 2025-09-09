#include "routines/ios.hpp"
#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "io/input.hpp"

#define EXTERN_MULTIGRID
# include "multiGridParams.hpp"

namespace nissa
{
  /// read the multigrid configuration file
  void read_DDalphaAMG_pars()
  {
    using namespace multiGrid;
    
    char path[]="DDalphaAMG_pars";
    
    smoother_iterations=4;
    
    nlevels=3;
    
    for(int ilev=0;ilev<nlevels;ilev++) nsetups[ilev]=4;
    for(int ilev=0;ilev<nlevels;ilev++) mu_factor[ilev]=1;
    for(int ilev=0;ilev<nlevels;ilev++) mu_factor_no_deflation[ilev]=1;
    for(int ilev=0;ilev<nlevels;ilev++) nu_pre[ilev]=0;
    for(int ilev=0;ilev<nlevels;ilev++) nu_post[ilev]=7;
    
    constexpr double def_coarse_solver_tol[MAX_MG_LEVELS]={0.15,0.22,0.46};
    for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
      {
	coarse_solver_tol_no_deflation[ilev]=
	  coarse_solver_tol[ilev]=
	  def_coarse_solver_tol[ilev];
      }
    constexpr int def_coarse_solver_maxiter[MAX_MG_LEVELS]={100,100,100};
    for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
      {
	coarse_solver_maxiter_no_deflation[ilev]=
	  coarse_solver_maxiter[ilev]=
	  def_coarse_solver_maxiter[ilev];
      }
    constexpr double def_smoother_tol[MAX_MG_LEVELS]={0.1,0.1,0.15};
    for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
      {
	smoother_tol_no_deflation[ilev]=
	  smoother_tol[ilev]=
	  def_smoother_tol[ilev];
      }
    constexpr double def_omega[MAX_MG_LEVELS]={0.85,0.85,0.85};
    for(int ilev=0;ilev<MAX_MG_LEVELS;ilev++)
      {
	omega_no_deflation[ilev]=
	  omega[ilev]=
	  def_omega[ilev];
      }
    
    
    if(nissa::fileExists(path))
      {
	nissa::open_input(path);
	int nr;
	
	do
	  {
	    char tag[100];
	    nr=nissa::read_var_catcherr(tag,"%s",100);
	    
	    int& nlevels=nissa::multiGrid::nlevels;
	    
	    if(nr>=1)
	      {
#define READ_VAR(TYPE,SHORT,NAME)					\
		if(strcasecmp(tag,#NAME)==0)				\
		  {							\
		    nissa::read_ ## TYPE(&NAME);			\
		    MASTER_PRINTF("DD: read " #NAME "=" SHORT "\n",NAME); \
		  }
		
#define READ_ARR(TYPE,SHORT,NAME)					\
		if(strcasecmp(tag,#NAME)==0)				\
		  for(int ilev=0;ilev<nlevels;ilev++)			\
		    {							\
		      nissa::read_ ## TYPE(&NAME[ilev]);		\
		      MASTER_PRINTF("DD: read " #NAME "[%d]=" SHORT "\n",ilev,NAME[ilev]); \
		    }
		
		READ_VAR(int,"%d",nlevels);
		READ_VAR(int,"%d",smoother_iterations);
		READ_VAR(double,"%lg",max_mass_for_deflation);
		READ_VAR(double,"%lg",max_mass);
		
		READ_ARR(int,"%d",nsetups);
		READ_ARR(double,"%lg",mu_factor);
		READ_ARR(double,"%lg",mu_factor_no_deflation);
		
		//size of the blocks, nb ought to historic bugs, one needs to provide block size in the order XTZY
		if(strcasecmp(tag,"block_size")==0)
		  {
		    block_size_set=true;
		    for(int ilev=0;ilev<nlevels;ilev++)
		      for(int idir=0;idir<4;idir++)
			{
			  int jdir=nissa::scidacMapping[idir];
			  nissa::read_int(&block_size[ilev][jdir]);
			  MASTER_PRINTF("DD: block_size[%d][%d*]=%d\n",ilev,jdir,block_size[ilev][jdir]);
			}
		    }
#ifdef USE_QUDA
		READ_VAR(double,"%lg",setup_refresh_tol);
		READ_VAR(double,"%lg",reliable_delta);
		READ_VAR(double,"%lg",reliable_delta_refinement);
		READ_VAR(int,"%d",nEigenvectors);
		READ_VAR(int,"%d",gcrNkrylov);
		READ_VAR(double,"%lg",eig_min);
		READ_VAR(double,"%lg",eig_max);
		READ_VAR(double,"%lg",setup_refresh_tol);
		
		READ_ARR(double,"%lg",omega);
		READ_ARR(double,"%lg",omega_no_deflation);
		READ_ARR(double,"%lg",coarse_solver_tol);
		READ_ARR(double,"%lg",coarse_solver_tol_no_deflation);
		READ_ARR(int,"%d",coarse_solver_maxiter);
		READ_ARR(int,"%d",coarse_solver_maxiter_no_deflation);
		READ_ARR(double,"%lg",smoother_tol);
		READ_ARR(double,"%lg",smoother_tol_no_deflation);
		READ_ARR(int,"%d",nu_pre);
		READ_ARR(int,"%d",nu_pre_no_deflation);
		READ_ARR(int,"%d",nu_post);
		READ_ARR(int,"%d",nu_post_no_deflation);
#endif
	      }
	    else MASTER_PRINTF("Finished reading the file '%s'\n",path);
	  }
	while(nr==1);
	
	nissa::close_input();
      }
    else MASTER_PRINTF("No '%s' file present, using standard configuration\n",path);
  }
  
  #undef READ_ARR
  #undef READ_VAR
}
