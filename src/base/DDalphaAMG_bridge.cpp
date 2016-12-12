#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <DDalphaAMG.h>

#include "base/debug.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "geometry/geometry_lx.hpp"
#include "io/checksum.hpp"
#include "new_types/su3.hpp"
#include "routines/thread.hpp"
#include "routines/mpi_routines.hpp"


#define EXTERN_DD_BRIDGE
#include "DDalphaAMG_bridge.hpp"

namespace DD
{
  int nlvl=2;
  DDalphaAMG_status status;
  int setup_status=0;
  
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
	bool import_as_diff=(check_old[i]==check_cur[i]);
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
      }
    else master_printf("DD: No import needed\n");
  }
  
  //initialize the bridge with DDalphaAMG
  void initialize(double kappa,double cSW,double mu)
  {
    static int inited=0;
    
    static DDalphaAMG_init init_params;
    static DDalphaAMG_parameters params;
    
    //check if cSW agreeing, when inited already
    if(inited and cSW!=init_params.csw)
      {
	master_printf("DD: cSW changed from %lg to %lg, reinitializing\n",init_params.csw,cSW);
	finalize();
	inited=0;
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
	for(int mu=0;mu<NDIM;mu++)
	  {
	    init_params.block_lattice[mu]=
	      (((nissa::glb_size[1]/nissa::nrank_dir[1])%2==0)?
	       (((nissa::glb_size[1]/nissa::nrank_dir[1])%4==0)?4:2):
	       (((nissa::glb_size[1]/nissa::nrank_dir[1])%3==0)?3:1));
	    init_params.theta[mu]=0;
	  }
	init_params.bc=0;
	
	//set the number of levels
	init_params.number_of_levels=nlvl;
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
	
	//parse success
	if(status.success!=nlvl)
	  {
	    nissa::master_fprintf(stderr,"WARNING: %d level initialized instead of %d\n",status.success,nlvl);
	    printf("DD WARNING: parameter: mg_lvl is changed to %d\n\n",status.success);
	    nlvl=status.success;
	  }
      }
    
    //to recheck
    DDalphaAMG_get_parameters(&params);
    
    //overwrite pars in any case
    params.print=1;
    //set transposer
    params.conf_index_fct=conf_index_fct;
    params.vector_index_fct=vector_index_fct;
    
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
	setup_status=0;
      }
    
    //update pars
    DDalphaAMG_update_parameters(&params,&status);
  }
  
  //setup DD if needed
  void update_setup()
  {
        //full setup
    if(setup_status==0)
      {
	master_printf("DD: Starting a new setup\n");
	DDalphaAMG_setup(&status);
      }
    
    setup_status=1;
  }
  
  //finalize
  void finalize()
  {
    master_printf("DD: finalizing\n");
    
    DDalphaAMG_finalize();
  }
  
  //solve
  int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision,nissa::spincolor *in)
  {
    initialize(kappa,cSW,mu);
    import_gauge_conf(conf);
    update_setup();
    
    //else DDalphaAMG_update_setup(int iterations, DDalphaAMG_status *mg_status)
    DDalphaAMG_solve((double*)out,(double*)in,precision,&status);
    //DDalphaAMG_apply_operator((double*)out,(double*)in,&status);
    
    master_printf("DD: Solving time %.2f sec (%.1f %% on coarse grid)\n",status.time,100.0*(status.coarse_time/status.time));
    master_printf("DD: Total iterations on fine grid %d\n", status.iter_count);
    master_printf("DD: Total iterations on coarse grids %d\n", status.coarse_iter_count);
    if(!status.success) crash("the solver did not converge!\n");
    //mul_r(phi_new ,mg_scale, phi_new, N);
    
    return status.success;
  }
}
