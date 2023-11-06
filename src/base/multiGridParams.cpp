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
    for(int ilev=0;ilev<nlevels;ilev++) nu_pre[ilev]=0;
    for(int ilev=0;ilev<nlevels;ilev++) nu_post[ilev]=7;
    
    if(nissa::file_exists(path))
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
		    master_printf("DD: read " #NAME "=" SHORT "\n",NAME); \
		  }
		
#define READ_ARR(TYPE,SHORT,NAME)					\
		if(strcasecmp(tag,#NAME)==0)				\
		  for(int ilev=0;ilev<nlevels;ilev++)			\
		    {							\
		      nissa::read_ ## TYPE(&NAME[ilev]);		\
		      master_printf("DD: read " #NAME "[%d]=" SHORT "\n",ilev,NAME[ilev]); \
		    }
		
		READ_VAR(int,"%d",nlevels);
		READ_VAR(int,"%d",smoother_iterations);
		READ_VAR(double,"%lg",max_mass_for_deflation);
		READ_VAR(double,"%lg",max_mass);
		
		READ_ARR(int,"%d",nsetups);
		READ_ARR(double,"%lg",mu_factor);
		
		//size of the blocks
		if(strcasecmp(tag,"block_size")==0)
		  {
		    block_size_set=true;
		    for(int ilev=0;ilev<nlevels;ilev++)
		      for(int idir=0;idir<4;idir++)
			{
			  int jdir=nissa::scidac_mapping[idir];
			  nissa::read_int(&block_size[ilev][jdir]);
			  master_printf("DD: block_size[%d][%d*]=%d\n",ilev,jdir,block_size[ilev][jdir]);
			}
		    }
#ifdef USE_QUDA
		READ_VAR(int,"%d",nEigenvectors);
		READ_VAR(int,"%d",gcrNkrylov);
		READ_VAR(double,"%lg",eig_min);
		READ_VAR(double,"%lg",eig_max);
		READ_ARR(int,"%d",nu_pre);
		READ_ARR(int,"%d",nu_post);
#endif
	      }
	    else master_printf("Finished reading the file '%s'\n",path);
	  }
	while(nr==1);
	
	nissa::close_input();
      }
    else master_printf("No '%s' file present, using standard configuration\n",path);
  }
  
  #undef READ_ARR
  #undef READ_VAR
}
