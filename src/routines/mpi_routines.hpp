#ifndef _MPI_ROUTINES_HPP
#define _MPI_ROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>
#include <algorithm>
#include <concepts>

#include <expr/comp.hpp>
#include <geometry/geometry_lx.hpp>
#include <metaprogramming/concepts.hpp>

#ifndef EXTERN_MPI
# define EXTERN_MPI extern
# define INIT_MPI_TO(...)
#else
# define INIT_MPI_TO(ARGS...) ARGS
#endif

namespace nissa
{
  
  
  /////////////////////////////////////////////////////////////////
  
#define DEFINE_MPI_DATATYPE_OF(T,MPI_T)		\
  /*! MPI Datatype corresponding to T */	\
  inline MPI_Datatype _MPI_Datatype_of(T*)	\
  {						\
    return MPI_T;				\
  }
  
  DEFINE_MPI_DATATYPE_OF(int64_t,MPI_LONG)
  DEFINE_MPI_DATATYPE_OF(double,MPI_DOUBLE)
  
  /// Instantiates the correct datatype, given the type
  template <typename T>
  MPI_Datatype MPI_Datatype_of()
  {
    return _MPI_Datatype_of((T*)nullptr);
  }
  
  /// Provides the correct operation for summing
  template <typename T>
  struct _MPI_Op_dispatcher
  {
    static MPI_Op sum()
    {
      return MPI_SUM;
    }
  };
  
#define DEFINE_MPI_OP_DISPATCHER(TYPE,OP)	\
  /*! Provides the correct operation for TYPE */\
  template <>					\
  struct _MPI_Op_dispatcher<TYPE>		\
  {						\
    /*! Sum operation */			\
    static MPI_Op sum()			\
    {						\
      return OP;				\
    }						\
  }
  
  /// Gets the sum operation for the type T
  template <typename T>
  MPI_Op MPI_Op_sum_for_type()
  {
   return _MPI_Op_dispatcher<T>::sum();
  }
  
#undef DEFINE_MPI_OP_DISPATCHER
  
  size_t MPI_Get_count_size_t(MPI_Status &status);
  void coords_broadcast(coords_t& c);
  void init_MPI_thread(int narg,char **arg);
  void ranks_abort(int err);
  void ranks_barrier();
  int broadcast(int in,int rank_from=0);
  double broadcast(double in,int rank_from=0);
#ifdef USE_MPI
  void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type);
#endif
  
  
  /////////////////////////////////////////////////////////////////
  
}

#endif
