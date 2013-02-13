#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

//to be moved elsewhere soon

#include <mpi.h>

//take the different with following multiple of eight
MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos)
{
  MPI_Offset diff=pos%8;
  if(diff!=0) diff=8-diff;

  return diff;
}

//ceil to next multiple of eight
MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos)
{return pos+diff_with_next_eight_multiple(pos);}

//summ two float128
void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type)
{for(int i=0;i<(*len);i++) float_128_summassign(((float_128*)out)[i],((float_128*)in)[i]);}

//broadcast an int
int master_broadcast(int in)
{
  MPI_Bcast(&in,1,MPI_INT,0,MPI_COMM_WORLD);
  return in;
}

//reduce a complex
void glb_reduce_complex(complex out_glb,complex in_loc)
{
#pragma omp single
  MPI_Allreduce(in_loc,reduce_complex,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  complex_copy(out_glb,reduce_complex);
}

//reduce a float_128
void glb_reduce_float_128(float_128 out_glb,float_128 in_loc)
{
#pragma omp single
  MPI_Allreduce(in_loc,reduce_float_128,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
  
  float_128_copy(out_glb,reduce_float_128);
}

//reduce a double
double glb_reduce_double(double in_loc)
{
#pragma omp single
  MPI_Allreduce(&in_loc,&reduce_double,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return reduce_double;
}

//reduce an int
int glb_reduce_int(int in_loc)
{
#pragma omp single
  MPI_Allreduce(&in_loc,&reduce_int,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  return reduce_int;
}

//reduce a double vector
void glb_reduce_double_vect(double *out_glb,double *in_loc,int nel)
{
#pragma omp single
  MPI_Allreduce(in_loc,out_glb,nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

//reduce a complex vector
void glb_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel)
{
#pragma omp single
  MPI_Allreduce((double*)in_loc,(double*)out_glb,2*nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}
