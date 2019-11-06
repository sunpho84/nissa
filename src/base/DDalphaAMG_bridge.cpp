#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <DDalphaAMG.h>

#include "base/debug.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "geometry/geometry_lx.hpp"
#include "io/checksum.hpp"
#include "io/input.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"


#define EXTERN_DD_BRIDGE
#include "DDalphaAMG_bridge.hpp"

namespace DD
{
  DDalphaAMG_init init_params;
  int &nlevels=init_params.number_of_levels;
  DDalphaAMG_parameters params;
  int nsetups[MAX_MG_LEVELS];
  int smoother_iterations;
  double mu_factor[MAX_MG_LEVELS];
  DDalphaAMG_status status;
  bool setup_valid=false;
  bool inited=false;
  bool block_size_set=false;
  nissa::coords block_size[MAX_MG_LEVELS];
  
  //remap swapping x and z
  void remap_coord(nissa::coords out,const nissa::coords in)
  {for(int mu=0;mu<NDIM;mu++) out[mu]=in[nissa::scidac_mapping[mu]];}
  
  //return the coordinate transposed
  static int cart_coords(MPI_Comm comm,int ext_rank,int maxdims,int c[])
  {
    nissa::coords int_c;
    int stat=MPI_Cart_coords(comm,ext_rank,maxdims,int_c);
    remap_coord(c,int_c);
    return stat;
  }
  
  //return the rank of remapped coordinate
  static int cart_rank(MPI_Comm comm,const int ext_c[],int *ext_rank)
  {
    nissa::coords c;
    remap_coord(c,ext_c);
    return MPI_Cart_rank(comm,c,ext_rank);
  }
  
  //return the index of the configuration
  static int conf_index_fct(int t,int z,int y,int x,int mu)
  {return sizeof(nissa::su3)/sizeof(double)*(nissa::scidac_mapping[mu]+NDIM*nissa::loclx_of_coord_list(t,x,y,z));}
  
  //return the index inside a spincolor
  static int vector_index_fct(int t,int z,int y,int x)
  {return sizeof(nissa::spincolor)/sizeof(double)*nissa::loclx_of_coord_list(t,x,y,z);}
  
  //read the nissa configuration file
  void read_DDalphaAMG_pars()
  {
    char path[]="DDalphaAMG_pars";
    
    smoother_iterations=4;
    
    nlevels=3;
    
    for(int ilev=0;ilev<nlevels;ilev++) nsetups[ilev]=4;
    for(int ilev=0;ilev<nlevels;ilev++) mu_factor[ilev]=1;
    
    if(nissa::file_exists(path))
      {
	nissa::open_input(path);
	int nr;
	
	do
	  {
	    char tag[100];
	    nr=nissa::read_var_catcherr(tag,"%s",100);
	    
	    if(nr>=1)
	      {
		//number of levels
		if(strcasecmp(tag,"nlevels")==0)
		  {
		    nissa::read_int(&nlevels);
		    master_printf("DD: read nlevels=%d\n",nlevels);
		  }
		//number of smoother iterations
		if(strcasecmp(tag,"smoother_iterations")==0)
		  {
		    nissa::read_int(&smoother_iterations);
		    master_printf("DD: read smoother_iterations=%d\n",smoother_iterations);
		  }
		//maximal mass
		if(strcasecmp(tag,"max_mass")==0)
		  {
		    nissa::read_double(&max_mass);
		    master_printf("DD: read max_mass=%lg\n",max_mass);
		  }
		//number of setups
		if(strcasecmp(tag,"nsetups")==0)
		  for(int ilev=0;ilev<nlevels;ilev++)
		    {
		      nissa::read_int(&nsetups[ilev]);
		      master_printf("DD: read nsetups[%d]=%d\n",ilev,nsetups[ilev]);
		    }
		//factor to increase mass in setup
		if(strcasecmp(tag,"mu_factor")==0)
		  for(int ilev=0;ilev<nlevels;ilev++)
		    {
		      nissa::read_double(&mu_factor[ilev]);
		      master_printf("DD: read mu_factor[%d]=%lg\n",ilev,mu_factor[ilev]);
		    }
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
	      }
	    else master_printf("Finished reading the file '%s'\n",path);
	  }
	while(nr==1);
	
	nissa::close_input();
      }
    else master_printf("No '%s' file present, using standard configuration\n",path);
  }
  
  //import a gauge configuration
  void import_gauge_conf(nissa::quad_su3 *conf)
  {
    static nissa::checksum check_old={0,0},check_cur;
    //compute checksum
    nissa::checksum_compute_nissa_data(check_cur,conf,sizeof(nissa::quad_su3),sizeof(double)*8);
    
    //verify if import needed
    bool import=false;
    for(int i=0;i<2;i++)
      {
	//check inited
	bool import_as_new=(check_old[i]==0);
	if(!import and import_as_new) master_printf("DD: Old checksum 0, need to import the conf\n");
	import|=import_as_new;
	//check diff
	bool import_as_diff=(check_old[i]!=check_cur[i]);
	if(!import and import_as_diff) master_printf("DD: Old checksum %d is %x, new is %x, need to import\n",i,check_old[i],check_cur[i]);
	import|=import_as_diff;
	//save
	check_old[i]=check_cur[i];
      }
    
    if(import)
      {
	DDalphaAMG_set_configuration((double*)conf,&status);
	if(status.success) verbosity_lv1_master_printf("DD: conf set, plaquette %e\n",status.info);
	else crash("configuration updating did not run correctly");
	setup_valid=false;
      }
    else master_printf("DD: No import needed\n");
  }
  
  //initialize the bridge with DDalphaAMG
  void initialize(double kappa,double cSW,double mu)
  {
    //check if cSW agreeing, when inited already
    if(inited and cSW!=init_params.csw)
      {
	master_printf("DD: cSW changed from %lg to %lg, reinitializing\n",init_params.csw,cSW);
	finalize();
      }
    
    if(!inited)
      {
	master_printf("DD: Not initialized, initializing\n");
	
	//communicator and transposers
	init_params.comm_cart=nissa::cart_comm;
	init_params.Cart_coords=cart_coords;
	init_params.Cart_rank=cart_rank;
	
	//sizes and coord
	remap_coord(init_params.global_lattice,nissa::glb_size);
	remap_coord(init_params.procs,nissa::nrank_dir);
	
	//block size and theta
	for(int dir=0;dir<NDIM;dir++)
	  {
	    int jdir=nissa::scidac_mapping[dir];
	    init_params.block_lattice[dir]=
	      (((nissa::glb_size[jdir]/nissa::nrank_dir[jdir])%2==0)?
	       (((nissa::glb_size[jdir]/nissa::nrank_dir[jdir])%4==0)?4:2):
	       (((nissa::glb_size[jdir]/nissa::nrank_dir[jdir])%3==0)?3:1));
	    if(block_size_set) init_params.block_lattice[dir]=block_size[0][dir];
	    master_printf("Dir %d block size: %d\n",dir,init_params.block_lattice[dir]);
	    init_params.theta[dir]=0;
	  }
	init_params.bc=0;
	
	//set threads
	init_params.number_openmp_threads=nissa::nthreads;
	//values for kappa, mu and csw
	init_params.kappa=kappa;
	init_params.mu=mu;
	init_params.csw=cSW;
	init_params.rnd_seeds=NULL;
	
	//init file
	init_params.init_file=NULL;
	
	//initialization
	DDalphaAMG_initialize(&init_params,&params,&status);
	inited=true;
      }
    
    //to recheck
    DDalphaAMG_get_parameters(&params);
    
    //block_size
    if(block_size_set)
      {
	block_size_set=false;
	
	for(int ilev=0;ilev<nlevels;ilev++)
	  for(int idir=0;idir<4;idir++)
	    params.block_lattice[ilev][idir]=block_size[ilev][idir];
      }
    
    //overwrite pars in any case
    params.print=1;
    //set transposer
    params.conf_index_fct=conf_index_fct;
    params.vector_index_fct=vector_index_fct;
    
    //check smoother_iterations
    master_printf("DD: Smoother iterations changed from %d to %d\n",params.smoother_iterations,smoother_iterations);
    params.smoother_iterations=smoother_iterations;
    
    //check mu_factor
    for(int ilev=0;ilev<nlevels;ilev++)
      if(mu_factor[ilev]!=params.mu_factor[ilev])
	{
	  master_printf("DD: Mu_factor for lev %d changed from %lg to %lg\n",ilev,params.mu_factor[ilev],mu_factor[ilev]);
	  params.mu_factor[ilev]=mu_factor[ilev];
	}
    
    //check nsetups
    for(int ilev=0;ilev<nlevels;ilev++)
      if(nsetups[ilev]!=params.setup_iterations[ilev])
	{
	  master_printf("DD: nsetups for lev %d changed from %d to %d\n",ilev,params.setup_iterations[ilev],nsetups[ilev]);
	  params.setup_iterations[ilev]=nsetups[ilev];
	}
    
    //check mass
    if(params.mu!=mu)
      {
	master_printf("DD: Mu changed from %lg to %lg\n",params.mu,mu);
	init_params.mu=params.mu=mu;
      }
    
    //check kappa
    if(params.kappa!=kappa)
      {
	master_printf("DD: kappa changed from %lg to %lg\n",params.kappa,kappa);
	init_params.kappa=params.kappa=kappa;
	setup_valid=false;
      }
    
    //update pars
    DDalphaAMG_update_parameters(&params,&status);
  }
  
  //setup DD if needed
  void update_setup()
  {
    //full setup
    if(!setup_valid)
      {
	master_printf("DD: Starting a new setup\n");
	DDalphaAMG_setup(&status);
      }
    
    setup_valid=true;
  }
  
  //finalize
  void finalize()
  {
    master_printf("DD: finalizing\n");
    
    DDalphaAMG_finalize();
    inited=false;
  }
  
  //solve
  int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in)
  {
    initialize(kappa,cSW,mu);
    import_gauge_conf(conf);
    update_setup();
    
    //else DDalphaAMG_update_setup(int iterations, DDalphaAMG_status *mg_status)
    DDalphaAMG_solve((double*)out,(double*)in,sqrt(precision2),&status);
    //DDalphaAMG_apply_operator((double*)out,(double*)in,&status);
    
    master_printf("DD: Solving time %.2f sec (%.1f %% on coarse grid)\n",status.time,100.0*(status.coarse_time/status.time));
    master_printf("DD: Total iterations on fine grid %d\n", status.iter_count);
    master_printf("DD: Total iterations on coarse grids %d\n", status.coarse_iter_count);
    if(!status.success) crash("the solver did not converge!\n");
    //mul_r(phi_new ,mg_scale, phi_new, N);
    
    return status.success;
  }
}
