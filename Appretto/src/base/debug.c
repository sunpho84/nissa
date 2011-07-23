#pragma once

#define fprintf_debug(fil,mes)						\
  MPI_Barrier(MPI_COMM_WORLD);						\
  if(rank==0) fprintf(fil,"Debug message: %s on file %s, line %d\n",mes,__FILE__,__LINE__);

#define printf_debug(mes) fprintf(stderr,mes);

#define crash(mes)				\
  {						\
    fprintf_debug(stderr,mes);			\
    if(rank==0) MPI_Abort(MPI_COMM_WORLD,1);	\
  }

#define free_null(a)							\
  {									\
    free(a);								\
    a=NULL;								\
  }

#define check_free(a)							\
  {									\
    if(a!=NULL) {free_null(a);}						\
    else crash("Error, trying to free an already freed array");		\
  }
