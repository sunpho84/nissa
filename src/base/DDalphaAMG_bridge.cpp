#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <DDalphaAMG.h>

#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/thread.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  namespace DD
  {
    int nlvl=2;
    DDalphaAMG_parameters params;
    DDalphaAMG_status status;
    
    //remap swapping x and z
    void remap_coord(coords out,const coords in)
    {for(int mu=0;mu<NDIM;mu++) out[mu]=in[map_mu[mu]];}
    
    //return the coordinate transposed
    static int cart_coords(MPI_Comm comm,int rank,int maxdims,int c[])
    {
      coords int_c;
      int stat=MPI_Cart_coords(comm,rank,maxdims,int_c);
      remap_coord(c,int_c);
      return stat;
    }
    
    //return the rank of remapped coordinate
    static int cart_rank(MPI_Comm comm,const int ext_c[],int *rank)
    {
      coords c;
      remap_coord(c,ext_c);
      return MPI_Cart_rank(comm,c,rank);
    }
    
    //return the index of the configuration
    static int conf_index_fct(int t,int z,int y,int x,int mu)
    {return sizeof(su3)/sizeof(double)*(map_mu[mu]+NDIM*loclx_of_coord_list(t,x,y,z));}
    
    //return the index inside a spincolor
    static int vector_index_fct(int t, int z, int y, int x )
    {return sizeof(spincolor)/sizeof(double)*loclx_of_coord_list(t,x,y,z);}
  }
  
  //initialize the bridge with DDalphaAMG
  void init_DDalphaAMG()
  {
    DDalphaAMG_init init;
    
    //communicator and transposers
    init.comm_cart=cart_comm;
    init.Cart_coords=DD::cart_coords;
    init.Cart_rank=DD::cart_rank;
    
    //sizes and coord
    DD::remap_coord(init.global_lattice,glb_size);
    DD::remap_coord(init.procs,nrank_dir);
    
    //block size
    for(int mu=0;mu<NDIM;mu++)
      {
	init.block_lattice[mu]=(((glb_size[1]/nrank_dir[1])%2==0)?
				(((glb_size[1]/nrank_dir[1])%4==0)?4:2):
				(((glb_size[1]/nrank_dir[1])%3==0)?3:1));
	init.theta[mu]=0;
      }
    
    //set the number of levels
    init.number_of_levels=DD::nlvl;
    //set threads
    init.number_openmp_threads=nthreads;
    //trivial values for kappa,mu and csw
    init.kappa=0.125;
    init.mu=0.5;
    init.csw=0;
    init.rnd_seeds=NULL;
    
    //init file
    init.init_file=NULL;
    
    //initialization
    DDalphaAMG_initialize(&init,&DD::params,&DD::status);
    
    //parse success
    if(DD::status.success!=DD::nlvl)
      {
	master_fprintf(stderr,"WARNING: %d level initialized instead of %d\n",DD::status.success,DD::nlvl);
	printf("MG WARNING: parameter: mg_lvl is changed to %d\n\n",DD::status.success);
	DD::nlvl=DD::status.success;
      }
    
    //set transposer
    DD::params.conf_index_fct=DD::conf_index_fct;
    DD::params.vector_index_fct=DD::vector_index_fct;
  }
}
