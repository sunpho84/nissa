#include "routines/ios.hpp"
#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_MULTIGRID
# include "multiGridParams.hpp"

namespace nissa
{
  /// read the multigrid configuration file
  void read_DDalphaAMG_pars()
  {
//     using namespace multiGrid;
    
//     char path[]="DDalphaAMG_pars";
    
//     smoother_iterations=4;
    
//     nlevels=3;
    
//     for(int ilev=0;ilev<nlevels;ilev++) nsetups[ilev]=4;
//     for(int ilev=0;ilev<nlevels;ilev++) mu_factor[ilev]=1;
//     for(int ilev=0;ilev<nlevels;ilev++) nu_pre[ilev]=0;
//     for(int ilev=0;ilev<nlevels;ilev++) nu_post[ilev]=7;
    
//     if(nissa::file_exists(path))
//       {
// 	nissa::open_input(path);
// 	int nr;
	
// 	do
// 	  {
// 	    char tag[100];
// 	    nr=nissa::read_var_catcherr(tag,"%s",100);
	    
// 	    int& nlevels=nissa::multiGrid::nlevels;
	    
// 	    if(nr>=1)
// 	      {
// 		//number of levels
// 		if(strcasecmp(tag,"nlevels")==0)
// 		  {
// 		    nissa::read_int(&nlevels);
// 		    master_printf("DD: read nlevels=%d\n",nlevels);
// 		  }
// 		//number of smoother iterations
// 		if(strcasecmp(tag,"smoother_iterations")==0)
// 		  {
// 		    nissa::read_int(&smoother_iterations);
// 		    master_printf("DD: read smoother_iterations=%d\n",smoother_iterations);
// 		  }
// 		//maximal mass
// 		if(strcasecmp(tag,"max_mass")==0)
// 		  {
// 		    nissa::read_double(&max_mass);
// 		    master_printf("DD: read max_mass=%lg\n",max_mass);
// 		  }
// 		//maximal mass to be used for deflation
// 		if(strcasecmp(tag,"max_mass_for_deflation")==0)
// 		  {
// 		    nissa::read_double(&max_mass_for_deflation);
// 		    master_printf("DD: read max_mass_for_deflation=%lg\n",max_mass_for_deflation);
// 		  }
// 		//number of setups
// 		if(strcasecmp(tag,"nsetups")==0)
// 		  for(int ilev=0;ilev<nlevels;ilev++)
// 		    {
// 		      nissa::read_int(&nsetups[ilev]);
// 		      master_printf("DD: read nsetups[%d]=%d\n",ilev,nsetups[ilev]);
// 		    }
// 		//factor to increase mass in setup
// 		if(strcasecmp(tag,"mu_factor")==0)
// 		  for(int ilev=0;ilev<nlevels;ilev++)
// 		    {
// 		      nissa::read_double(&mu_factor[ilev]);
// 		      master_printf("DD: read mu_factor[%d]=%lg\n",ilev,mu_factor[ilev]);
// 		    }
// 		//size of the blocks
// 		if(strcasecmp(tag,"block_size")==0)
// 		  {
// 		    block_size_set=true;
// 		    for(int ilev=0;ilev<nlevels;ilev++)
// 		      for(int idir=0;idir<4;idir++)
// 			{
// 			  int jdir=nissa::scidac_mapping[idir];
// 			  nissa::read_int(&block_size[ilev][jdir]);
// 			  master_printf("DD: block_size[%d][%d*]=%d\n",ilev,jdir,block_size[ilev][jdir]);
// 			}
// 		    }
// #ifdef USE_QUDA
// 		//orthogonalization before
// 		if(strcasecmp(tag,"nu_pre")==0)
// 		  for(int ilev=0;ilev<nlevels;ilev++)
// 		    {
// 		      nissa::read_int(&nu_pre[ilev]);
// 		      master_printf("DD: read nu_pre[%d]=%d\n",ilev,nu_pre[ilev]);
// 		    }
// 		//orthogonalization before
// 		if(strcasecmp(tag,"nu_post")==0)
// 		  for(int ilev=0;ilev<nlevels;ilev++)
// 		    {
// 		      nissa::read_int(&nu_post[ilev]);
// 		      master_printf("DD: read nu_post[%d]=%d\n",ilev,nu_post[ilev]);
// 		    }
// 		//number of eigenvectors
// 		if(strcasecmp(tag,"nEigenvectors")==0)
// 		  {
// 		    nissa::read_int(&nEigenvectors);
// 		    master_printf("DD: nEigenvectors=%d\n",nEigenvectors);
// 		  }
// 		//minimal eigenvalue
// 		if(strcasecmp(tag,"eig_min")==0)
// 		  {
// 		    nissa::read_double(&eig_min);
// 		    master_printf("DD: read eig_min=%lg\n",eig_min);
// 		  }
// 		//maximal eigenvalue
// 		if(strcasecmp(tag,"eig_max")==0)
// 		  {
// 		    nissa::read_double(&eig_max);
// 		    master_printf("DD: read eig_max=%lg\n",eig_max);
// 		  }
// #endif
// 	      }
// 	    else master_printf("Finished reading the file '%s'\n",path);
// 	  }
// 	while(nr==1);
	
// 	nissa::close_input();
//       }
//     else master_printf("No '%s' file present, using standard configuration\n",path);
    
    crash("");
  }
}
