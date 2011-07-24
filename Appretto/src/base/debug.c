#pragma once

#define fprintf_debug(fil,mes)						\
  {									\
    MPI_Barrier(MPI_COMM_WORLD);					\
    if(rank==0)								\
      {									\
	fprintf(fil,"Debug message: %s on file %s, line %d\n",mes,__FILE__,__LINE__); \
	fflush(fil);							\
      }									\
    MPI_Barrier(MPI_COMM_WORLD);					\
  }

#define printf_debug(mes) fprintf_debug(stdout,mes);

#define crash(mes)				\
  {						\
    fprintf_debug(stderr,mes);			\
    if(rank==0) MPI_Abort(MPI_COMM_WORLD,1);	\
  }

void free_null(void *a)
{
  free(a);
  a=NULL;
}

void check_free(void *a)
  {
    if(a!=NULL) free_null(a);
    else crash("Error, trying to free an already freed array");
  }
