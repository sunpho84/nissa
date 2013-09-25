#ifndef _MPI_NISSA_H
#define _MPI_NISSA_H

#ifdef USE_MPI
 #include <mpi.h>
#endif

void coords_broadcast(coords c);
void get_MPI_nranks();
void get_MPI_rank();
void init_MPI_thread(int narg,char **arg);
void define_MPI_types();
void create_MPI_cartesian_grid();
void ranks_abort(int err);
void ranks_barrier();
double glb_reduce_double(double in_loc);
int glb_reduce_int(int in_loc);
int master_broadcast(int in);
#ifdef USE_MPI
MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos);
MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos);
void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type);
#else
uint64_t ceil_to_next_eight_multiple(uint64_t pos);
uint64_t diff_with_next_eight_multiple(uint64_t pos);
#endif
void glb_reduce_complex(complex out_glb,complex in_loc);
void glb_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel);
void glb_reduce_double_vect(double *out_glb,double *in_loc,int nel);
void glb_reduce_float_128(float_128 out_glb,float_128 in_loc);

#endif
