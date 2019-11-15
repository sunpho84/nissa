#ifndef _MPI_ROUTINES_HPP
#define _MPI_ROUTINES_HPP

#include <mpi.h>
#include <algorithm>

#include "geometry/geometry_lx.hpp"
#include "math_routines.hpp"
#include "new_types/float_128.hpp"
#include "new_types/rat_approx.hpp"

#ifndef EXTERN_MPI
 #define EXTERN_MPI extern
#endif

namespace nissa
{
  //basic mpi types
  EXTERN_MPI MPI_Datatype MPI_FLOAT_128;
  EXTERN_MPI MPI_Datatype MPI_SU3;
  EXTERN_MPI MPI_Datatype MPI_QUAD_SU3;
  EXTERN_MPI MPI_Datatype MPI_AS2T_SU3;
  EXTERN_MPI MPI_Datatype MPI_COLOR;
  EXTERN_MPI MPI_Datatype MPI_SPIN;
  EXTERN_MPI MPI_Datatype MPI_SPINSPIN;
  EXTERN_MPI MPI_Datatype MPI_SPINCOLOR;
  EXTERN_MPI MPI_Datatype MPI_SPINCOLOR_128;
  EXTERN_MPI MPI_Datatype MPI_REDSPINCOLOR;
  //float 128 summ
  EXTERN_MPI MPI_Op MPI_FLOAT_128_SUM;
  
  EXTERN_MPI MPI_Datatype MPI_LX_SU3_EDGES_SEND[NDIM*(NDIM-1)/2],MPI_LX_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_LX_AS2T_SU3_EDGES_SEND[NDIM*(NDIM-1)/2],MPI_LX_AS2T_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_LX_QUAD_SU3_EDGES_SEND[NDIM*(NDIM-1)/2],MPI_LX_QUAD_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_EO_QUAD_SU3_EDGES_SEND[96],MPI_EO_QUAD_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  
  //volume, plan and line communicator
  EXTERN_MPI MPI_Comm cart_comm;
  EXTERN_MPI MPI_Comm plan_comm[NDIM];
  EXTERN_MPI MPI_Comm line_comm[NDIM];
  
  size_t MPI_Get_count_size_t(MPI_Status &status);
  void coords_broadcast(coords c);
  void get_MPI_nranks();
  void get_MPI_rank();
  void init_MPI_thread(int narg,char **arg);
  void define_MPI_types();
  void create_MPI_cartesian_grid();
  void ranks_abort(int err);
  void ranks_barrier();
  // float glb_reduce_single(float in_loc);
  // double glb_reduce_double(double in_loc,double (*thread_op)(double,double)=summ<double>,MPI_Op mpi_op=MPI_SUM);
  // inline double glb_max_double(double in_loc)
  // {return glb_reduce_double(in_loc,nissa_max,MPI_MAX);}
  int broadcast(int in,int rank_from=0);
  void broadcast(rat_approx_t *rat,int rank_from=0);
  double broadcast(double in,int rank_from=0);
#ifdef USE_MPI
  MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos);
  MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos);
  void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type);
#else
  uint64_t ceil_to_next_eight_multiple(uint64_t pos);
  uint64_t diff_with_next_eight_multiple(uint64_t pos);
#endif
  
  /////////////////////////////////////////////////////////////////
  
  //reduce a double
  inline double MPI_reduce_double(double in_loc,MPI_Op mpi_op)
  {
    double out_glb;
    
    MPI_Allreduce(&in_loc,&out_glb,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce a double
  inline float MPI_reduce_single(float in_loc)
  {
    float out_glb;
    
    MPI_Allreduce(&in_loc,&out_glb,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce an int
  inline int MPI_reduce_int(int in_loc)
  {
    int out_glb;
    
    MPI_Allreduce(&in_loc,&out_glb,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce a complex
  inline void MPI_reduce_complex(complex out_glb,complex in_loc)
  {
    MPI_Allreduce(&in_loc,&out_glb,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }
  
  //reduce a float_128
  inline float_128 MPI_reduce_float_128(const float_128& in_loc)
  {
    float_128 out_glb;
    
    MPI_Allreduce(&in_loc,&out_glb,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce a complex 128
  inline void MPI_reduce_complex_128(complex_128 out_glb,complex_128 in_loc)
  {
    MPI_Allreduce(&in_loc,&out_glb,2,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
  }
  
  //reduce a double vector
  inline void MPI_reduce_double_vect(double *out_glb,double *in_loc,int nel)
  {
    MPI_Allreduce(in_loc,out_glb,nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }
  
  inline void MPI_reduce_double_vect(double *vect,int nel)
  {
    MPI_reduce_double_vect(vect,(double*)MPI_IN_PLACE,nel);
  }
  
  inline void MPI_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel)
  {
    MPI_reduce_double_vect(out_glb[0],in_loc[0],2*nel);
  }
  
  inline void MPI_reduce_complex_vect(complex *vect,int nel)
  {
    MPI_reduce_double_vect(vect[0],2*nel);
  }
}

#endif
