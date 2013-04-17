#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>

#include "../base/global_variables.h"
#include "../base/thread_macros.h"
#include "../new_types/complex.h"
#include "../new_types/float128.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/thread.h"

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

//reduce a double
double glb_reduce_double(double in_loc)
{
  double out_glb;
  
  if(thread_pool_locked) MPI_Allreduce(&in_loc,&out_glb,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  else
    {
      GET_THREAD_ID();
      
      //copy loc in the buf and sync all the threads
      glb_double_reduction_buf[thread_id]=in_loc;
      thread_barrier(DOUBLE_REDUCE_FIRST_BARRIER);
      
      //within master thread summ all the pieces and between MPI
      if(IS_MASTER_THREAD)
	{
	  for(int ith=1;ith<nthreads;ith++) in_loc+=glb_double_reduction_buf[ith];
	  MPI_Allreduce(&in_loc,&(glb_double_reduction_buf[0]),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#pragma omp flush
	}
      
      //read glb val
      THREAD_ATOMIC_EXEC(out_glb=glb_double_reduction_buf[0];);
    }
  
  return out_glb;
}

//reduce an int
void glb_reduce_int(int *out_glb,int in_loc)
{
  if(thread_pool_locked) MPI_Allreduce(&in_loc,out_glb,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  else crash("not threaded yet");
}

//reduce a complex
void glb_reduce_complex(complex out_glb,complex in_loc)
{for(int ri=0;ri<2;ri++) out_glb[ri]=glb_reduce_double(in_loc[ri]);}

//reduce a float_128
void glb_reduce_float_128(float_128 out_glb,float_128 in_loc)
{
  if(thread_pool_locked) MPI_Allreduce(in_loc,out_glb,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
  else
    {
      GET_THREAD_ID();
      
      //copy loc in the buf and sync all the threads
      float_128_copy(glb_float_128_reduction_buf[thread_id],in_loc);
      thread_barrier(FLOAT_128_REDUCE_FIRST_BARRIER);
      
      //within master thread summ all the pieces and between MPI
      if(IS_MASTER_THREAD)
	{
	  for(int ith=1;ith<nthreads;ith++) float_128_summassign(in_loc,glb_float_128_reduction_buf[ith]);
	  MPI_Allreduce(in_loc,glb_float_128_reduction_buf[0],1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
#pragma omp flush
	}
      
      //read glb val
      THREAD_ATOMIC_EXEC(float_128_copy(out_glb,glb_float_128_reduction_buf[0]););
    }
}

//reduce a complex 128
void glb_reduce_complex_128(complex_128 out_glb,complex_128 in_loc)
{for(int ri=0;ri<2;ri++) glb_reduce_float_128(out_glb[ri],in_loc[ri]);}

//reduce a double vector
void glb_reduce_double_vect(double *out_glb,double *in_loc,int nel)
{MPI_Allreduce(in_loc,out_glb,nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}

//reduce a complex vector
void glb_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel)
{MPI_Allreduce((double*)in_loc,(double*)out_glb,2*nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}
