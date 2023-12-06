#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

/// \file routines/mpi_routines.cpp

#ifdef USE_MPI
# include <mpi.h>
#endif

#define EXTERN_MPI
# include <routines/mpi_routines.hpp>

#include <geometry/geometry_lx.hpp>
#include <linalgs/reduce.hpp>
#include <routines/ios.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  //return  the count covnerted to size_t
  size_t MPI_Get_count_size_t(MPI_Status &status)
  {
    int nbytes;
    decript_MPI_error(MPI_Get_count(&status,MPI_BYTE,&nbytes),"while counting bytes");
    if(nbytes<0) crash("negative count: %d",nbytes);
    
    return (size_t)nbytes;
  }
  
  //take the different with following multiple of eight
  MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos)
  {
    MPI_Offset diff=pos%8;
    if(diff!=0) diff=8-diff;
    
    return diff;
  }
  
  //init mpi
  void init_MPI_thread(int narg,char **arg)
  {
#ifdef USE_MPI
    
# ifdef USE_OPENMP
    int provided;
    MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
# else
    MPI_Init(&narg,&arg);
 #endif
#endif
  }
  
  //define the cartesian grid
  void create_MPI_cartesian_grid()
  {
#ifdef USE_MPI
    rank_coord=coord_of_rank(thisRank());
#else
    for(int mu=0;mu<NDIM;mu++) rank_coord[mu]=0;
#endif
  }
  
  //barrier
  void ranks_barrier()
  {
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  //abort
  void ranks_abort(int err)
  {
#ifdef USE_MPI
    printf("on rank %ld aborting\n",thisRank());
    MPI_Abort(MPI_COMM_WORLD,0);
#else
    exit(0);
#endif
  }

  //broadcast a coord
  void coords_broadcast(coords_t& c)
  {
    MPI_Bcast(&c[0],NDIM,MPI_INT,thisRank(),MPI_COMM_WORLD);
  }
  
  //ceil to next multiple of eight
  MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos)
  {return pos+diff_with_next_eight_multiple(pos);}
  
  //internal version
  template <class T>
  T broadcast_internal(T in,int rank_from,MPI_Datatype type)
  {
    T out=in;
    MPI_Bcast(&out,1,type,rank_from,MPI_COMM_WORLD);
    
    return out;
  }
  
  //broadcast an int
  int broadcast(int in,int rank_from)
  {return broadcast_internal(in,rank_from,MPI_INT);}
  
  //broadcast an int
  double broadcast(double in,int rank_from)
  {return broadcast_internal(in,rank_from,MPI_DOUBLE);}
}
