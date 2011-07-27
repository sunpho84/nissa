#pragma once

#define debug 1

#define appretto_vect_string_length 20

#define appretto_malloc(a,b,c) appretto_true_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define appretto_free(a) appretto_true_free((void**)&(a),__FILE__,__LINE__)

#define fprintf_debug(fil,mes)                                          \
  {                                                                     \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    if(rank==0)                                                         \
      {                                                                 \
        fprintf(fil,"Debug message: %s on file %s, line %d\n",mes,__FILE__,__LINE__); \
        fflush(fil);                                                    \
      }                                                                 \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
  }

#define printf_debug(mes) fprintf_debug(stdout,mes);

#define crash(mes)                              \
  {                                             \
    fprintf_debug(stderr,mes);                  \
    if(rank==0) MPI_Abort(MPI_COMM_WORLD,1);    \
  }
