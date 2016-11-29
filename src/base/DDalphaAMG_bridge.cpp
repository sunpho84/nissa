#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <DDalphaAMG.h>

#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/thread.hpp"
#include "routines/mpi_routines.hpp"


#define EXTERN_DD_BRIDGE
#include "DDalphaAMG_bridge.hpp"

namespace DD
{
  int nlvl=2;
  DDalphaAMG_parameters params;
  DDalphaAMG_status status;
  
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
    for(int mu=0;mu<NDIM;mu++) c[mu]=(c[mu]+nissa::nrank_dir[mu])%nissa::nrank_dir[mu];
    return MPI_Cart_rank(comm,c,ext_rank);
  }
  
  //return the index of the configuration
  static int conf_index_fct(int t,int z,int y,int x,int mu)
  {return sizeof(nissa::su3)/sizeof(double)*(nissa::scidac_mapping[mu]+NDIM*nissa::loclx_of_coord_list(t,x,y,z));}
  
  //return the index inside a spincolor
  static int vector_index_fct(int t, int z, int y, int x )
  {return sizeof(nissa::spincolor)/sizeof(double)*nissa::loclx_of_coord_list(t,x,y,z);}
  
  //import a gauge configuration
  void import_gauge_conf(nissa::quad_su3 *conf)
  {
    DDalphaAMG_set_configuration((double*)conf,&status);
    if(status.success) verbosity_lv1_master_printf("DDalphaAMG conf set, plaquette %e\n",status.info);
    else crash("configuration updating did not run correctly");
  }
  
  //set the value of the mass
  void set_mass(double mass)
  {
    DDalphaAMG_get_parameters(&params);
    if(params.mu!=mass)
      {
	params.mu=mass;
	DDalphaAMG_update_parameters(&params,&status);
      }
  }
  
  //initialize the bridge with DDalphaAMG
  void init_DDalphaAMG()
  {
    DDalphaAMG_init init;
    
    //communicator and transposers
    init.comm_cart=nissa::cart_comm;
    init.Cart_coords=cart_coords;
    init.Cart_rank=cart_rank;
    
    //sizes and coord
    remap_coord(init.global_lattice,nissa::glb_size);
    remap_coord(init.procs,nissa::nrank_dir);
    
    //block size and theta
    for(int mu=0;mu<NDIM;mu++)
      {
	init.block_lattice[mu]=
	  (((nissa::glb_size[1]/nissa::nrank_dir[1])%2==0)?
	   (((nissa::glb_size[1]/nissa::nrank_dir[1])%4==0)?4:2):
	   (((nissa::glb_size[1]/nissa::nrank_dir[1])%3==0)?3:1));
	init.theta[mu]=0;
      }
    init.bc=0;
    
    //set the number of levels
    init.number_of_levels=nlvl;
    //set threads
    init.number_openmp_threads=nissa::nthreads;
    //trivial values for kappa,mu and csw
    init.kappa=0.125;
    init.mu=0.5;
    init.csw=0;
    init.rnd_seeds=NULL;
    
    //init file
    init.init_file=NULL;
    
    //initialization
    DDalphaAMG_initialize(&init,&params,&status);
    
    //parse success
    if(status.success!=nlvl)
      {
	nissa::master_fprintf(stderr,"WARNING: %d level initialized instead of %d\n",status.success,nlvl);
	printf("MG WARNING: parameter: mg_lvl is changed to %d\n\n",status.success);
	nlvl=status.success;
      }
    
    //set transposer
    params.print=1;
    params.conf_index_fct=conf_index_fct;
    params.vector_index_fct=vector_index_fct;
    DDalphaAMG_update_parameters(&params,&status);
  }
  
  //finalize
  void finalize_DDalphaAMG()
  {DDalphaAMG_finalize();}
  
  //solve
  int solve(nissa::spincolor *out,nissa::spincolor *in,double mu,double precision)
  {
    set_mass(mu);
    static int setup=1;
    if(setup)
      {
	DDalphaAMG_setup(&status);
	setup=0;
      }
    //else DDalphaAMG_update_setup(int iterations, DDalphaAMG_status *mg_status)
    DDalphaAMG_solve((double*)out,(double*)in,precision,&status);
    
    master_printf("Solving time %.2f sec (%.1f %% on coarse grid)\n",status.time,100.0*(status.coarse_time/status.time));
    master_printf("Total iterations on fine grid %d\n", status.iter_count);
    master_printf("Total iterations on coarse grids %d\n", status.coarse_iter_count);
    if(!status.success) crash("the solver did not converge!\n");
    //mul_r(phi_new ,mg_scale, phi_new, N);
    
    return status.success;
  }
}
