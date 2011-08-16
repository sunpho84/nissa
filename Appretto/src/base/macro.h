#pragma once

#define debug 1

#define appretto_vect_string_length 20

#define appretto_malloc(a,b,c) appretto_true_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define appretto_free(a) a=appretto_true_free(a,__FILE__,__LINE__)

#define fprintf_debug(fil,...)					\
  {                                                                     \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    if(rank==0)                                                         \
      {                                                                 \
        fprintf(fil,"Debug message: <<");				\
	fprintf(fil,__VA_ARGS__);					\
	fprintf(fil,">> on file %s, line %d\n",__FILE__,__LINE__);	\
	fflush(fil);							\
      }                                                                 \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
  }

#define printf_debug(...) fprintf_debug(stdout,__VA_ARGS__);

#define crash(...)                              \
  {                                             \
    fprintf_debug(stderr,__VA_ARGS__);		\
    fflush(stderr);				\
    if(rank==0) MPI_Abort(MPI_COMM_WORLD,1);    \
  }

#define if_appretto_vect_not_initialized() if(main_appretto_arr!=((char*)&main_appretto_vect)+sizeof(appretto_vect))

