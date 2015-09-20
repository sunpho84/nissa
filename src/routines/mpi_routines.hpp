#ifndef _MPI_NISSA_HPP
#define _MPI_NISSA_HPP

#include <mpi.h>
#include <algorithm>
#include "new_types/new_types_definitions.hpp"
#include "math_routines.hpp"

namespace nissa
{
  size_t MPI_Get_count_size_t(MPI_Status &status);
  void coords_broadcast(coords c);
  void get_MPI_nranks();
  void get_MPI_rank();
  void init_MPI_thread(int narg,char **arg);
  void define_MPI_types();
  void create_MPI_cartesian_grid();
  void ranks_abort(int err);
  void ranks_barrier();
  float glb_reduce_single(float in_loc);
  double glb_reduce_double(double in_loc,double (*thread_op)(double,double)=summ<double>,MPI_Op mpi_op=MPI_SUM);
  inline double glb_max_double(double in_loc)
  {return glb_reduce_double(in_loc,nissa_max,MPI_MAX);}
  int glb_reduce_int(int in_loc);
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
  void glb_reduce_complex(complex out_glb,complex in_loc);

  void glb_nodes_reduce_double_vect(double *out_glb,double *in_loc,int nel);
  inline void glb_nodes_reduce_double_vect(double *vect,int nel)
  {glb_nodes_reduce_double_vect(vect,(double*)MPI_IN_PLACE,nel);}
  inline void glb_nodes_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel)
  {glb_nodes_reduce_double_vect(out_glb[0],in_loc[0],2*nel);}
  inline void glb_nodes_reduce_complex_vect(complex *vect,int nel)
  {glb_nodes_reduce_double_vect(vect[0],2*nel);}
  void glb_reduce_float_128(float_128 out_glb,float_128 in_loc);
}
#endif
