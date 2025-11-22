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
		
		using namespace nissa::multiGrid::internal;
		
#define READ_ARR_OPT_DEFL(TYPE,SHORT,NAME)				\
		READ_ARR(TYPE,SHORT,NAME);				\
		READ_ARR(TYPE,SHORT,NAME ## _no_deflation)
		
#define READ_VAR_OPT_DEFL(TYPE,SHORT,NAME)				\
		READ_VAR(TYPE,SHORT,NAME);				\
		READ_VAR(TYPE,SHORT,NAME ## _no_deflation)
		
		READ_VAR(int,"%d",nlevels);
		READ_VAR(int,"%d",smoother_iterations);
		READ_VAR(double,"%lg",max_mass_for_deflation);
		READ_VAR(double,"%lg",max_mass);
		
		READ_ARR(int,"%d",nsetups);
		
		READ_ARR_OPT_DEFL(double,"%lg",mu_factor);
		
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
		READ_VAR(int,"%d",nEigenvectors);
		READ_VAR(int,"%d",gcrNkrylov);
		READ_VAR(double,"%lg",eig_min);
		READ_VAR(double,"%lg",eig_max);
		READ_VAR(double,"%lg",setup_refresh_tol);
		
		READ_VAR_OPT_DEFL(double,"%lg",reliable_delta);
		READ_VAR_OPT_DEFL(double,"%lg",reliable_delta_refinement);
		
		READ_ARR_OPT_DEFL(int,"%d",nu_post);
		READ_ARR_OPT_DEFL(double,"%lg",coarse_solver_tol);
		READ_ARR_OPT_DEFL(int,"%d",coarse_solver_maxiter);
		READ_ARR_OPT_DEFL(double,"%lg",smoother_tol);
		READ_ARR_OPT_DEFL(int,"%d",nu_pre);
		READ_ARR_OPT_DEFL(int,"%d",nu_post);
		READ_ARR_OPT_DEFL(double,"%lg ",omega);
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
  #undef READ_ARR_OPT_DEFL
  #undef READ_VAR_OPT_DEFL
}
