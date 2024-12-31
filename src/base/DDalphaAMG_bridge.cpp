#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_DD_BRIDGE
# include "DDalphaAMG_bridge.hpp"

#include "base/debug.hpp"
#include "base/export_conf_to_external_solver.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "geometry/geometry_lx.hpp"
#include "io/checksum.hpp"
#include "io/input.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

namespace DD
{
  DDalphaAMG_init init_params;
  DDalphaAMG_parameters params;
  bool inited=false;
  
  //remap swapping x and z
  void remap_coord(int* out,const int* in)
  {for(int mu=0;mu<NDIM;mu++) out[mu]=in[nissa::scidacMapping[mu]];}
  
  //return the coordinate transposed
  static int cart_coords(MPI_Comm comm,int ext_rank,int maxdims,int c[])
  {
    int int_c[NDIM];
    int stat=MPI_Cart_coords(comm,ext_rank,maxdims,int_c);
    remap_coord(c,int_c);
    return stat;
  }
  
  //return the rank of remapped coordinate
  static int cart_rank(MPI_Comm comm,const int ext_c[],int *ext_rank)
  {
    int c[NDIM];
    remap_coord(c,ext_c);
    return MPI_Cart_rank(comm,c,ext_rank);
  }
  
  //return the index of the configuration
  static int conf_index_fct(int t,int z,int y,int x,int mu)
  {return sizeof(nissa::su3)/sizeof(double)*(nissa::scidacMapping[mu]+NDIM*nissa::loclxOfCoordList(t,x,y,z));}
  
  //return the index inside a spincolor
  static int vector_index_fct(int t,int z,int y,int x)
  {return sizeof(nissa::spincolor)/sizeof(double)*nissa::loclxOfCoordList(t,x,y,z);}
  
  /// Check if cSW is changed
  bool check_cSW_changed(const double& cSW)
  {
    /// Check if cSW is changed
    const bool cSW_changed=(inited and cSW!=init_params.csw);
    
    if(cSW_changed)
      master_printf("DD: cSW changed from %lg to %lg\n",init_params.csw,cSW);
    
    return cSW_changed;
  }
  
  /// Check if kappa or cSW are changed
  bool check_kappa_changed(const double& kappa)
  {
    /// Check if kappa is changed
    const bool kappa_changed=(inited and kappa!=init_params.kappa);
    
    if(kappa_changed)
      master_printf("DD: kappa changed from %lg to %lg\n",init_params.kappa,kappa);
    
    return kappa_changed;
  }
  
  /// Initialize the bridge with DDalphaAMG, for the given kappa and cSW
  void initialize(const double& kappa,const double& cSW,const double& mu)
  {
    if(check_kappa_changed(kappa) or check_cSW_changed(cSW))
      {
	master_printf("Kappa or cSW changed, reinitializing\n");
	finalize();
      }
    
    if(not inited)
      {
	master_printf("DD: Not initialized, initializing\n");
	
	//communicator and transposers
	crash("reimplement");
	// init_params.comm_cart=nissa::cart_comm;
	init_params.Cart_coords=cart_coords;
	init_params.Cart_rank=cart_rank;
	init_params.number_of_levels=nissa::multiGrid::nlevels;
	
	//sizes and coord
	remap_coord(init_params.global_lattice,&nissa::glbSize[0]);
	remap_coord(init_params.procs,&nissa::nRanksDir[0]);
	
	//block size and theta
	for(int dir=0;dir<NDIM;dir++)
	  {
	    int jdir=nissa::scidacMapping[dir];
	    init_params.block_lattice[dir]=
	      (((nissa::glbSize[jdir]/nissa::nRanksDir[jdir])%2==0)?
	       (((nissa::glbSize[jdir]/nissa::nRanksDir[jdir])%4==0)?4:2):
	       (((nissa::glbSize[jdir]/nissa::nRanksDir[jdir])%3==0)?3:1));
	    if(nissa::multiGrid::block_size_set)
	      init_params.block_lattice[dir]=nissa::multiGrid::block_size[0][dir];
	    master_printf("Dir %d block size: %d\n",dir,init_params.block_lattice[dir]);
	    init_params.theta[dir]=0;
	  }
	init_params.bc=0;
	
	//set threads
	//init_params.number_openmp_threads=nissa::nthreads;
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
    
    using namespace nissa::multiGrid;
    
    //block_size
    if(block_size_set)
      for(int ilev=0;ilev<nlevels;ilev++)
	for(int idir=0;idir<4;idir++)
	  params.block_lattice[ilev][idir]=block_size[ilev][idir];
    
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
  
  void set_configuration(nissa::quad_su3* conf)
  {
    DDalphaAMG_set_configuration((double*)conf,&DD::status);
  }
  
  //setup DD if needed
  void update_setup()
  {
    bool& setup_valid=nissa::multiGrid::setup_valid;
    
    //full setup
    if(not setup_valid)
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
  
  //solve: if squared is asked, solve the odd square (Dkern_squared)
  int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in,const bool squared)
  {
    initialize(kappa,cSW,mu);
    nissa::export_gauge_conf_to_external_solver(conf);
    update_setup();
    
    //else DDalphaAMG_update_setup(int iterations, DDalphaAMG_status *mg_status)
    if(squared) DDalphaAMG_solve_squared_odd((double*)out,(double*)in,sqrt(precision2),&status);
    else        DDalphaAMG_solve((double*)out,(double*)in,sqrt(precision2),&status);
    //DDalphaAMG_apply_operator((double*)out,(double*)in,&status);
    
    master_printf("DD: Solving time %.2f sec (%.1f %% on coarse grid)\n",status.time,100.0*(status.coarse_time/status.time));
    master_printf("DD: Total iterations on fine grid %d\n", status.iter_count);
    master_printf("DD: Total iterations on coarse grids %d\n", status.coarse_iter_count);
    if(!status.success) crash("the solver did not converge!\n");
    //mul_r(phi_new ,mg_scale, phi_new, N);
    
    return status.success;
  }
}
