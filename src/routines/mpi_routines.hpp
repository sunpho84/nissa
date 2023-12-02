#ifndef _MPI_ROUTINES_HPP
#define _MPI_ROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>
#include <algorithm>
#include <concepts>

#include "expr/baseComp.hpp"
#include "geometry/geometry_lx.hpp"
#include "math_routines.hpp"
#include "metaprogramming/concepts.hpp"

#ifndef EXTERN_MPI
# define EXTERN_MPI extern
# define INIT_MPI_TO(...)
#else
# define INIT_MPI_TO(ARGS...) ARGS
#endif

namespace nissa
{
  
  EXTERN_MPI int master_rank INIT_MPI_TO(=0);
  
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
  void get_MPI_nranks();
  void get_MPI_rank();
  void init_MPI_thread(int narg,char **arg);
  void create_MPI_cartesian_grid();
  void ranks_abort(int err);
  void ranks_barrier();
  int broadcast(int in,int rank_from=0);
  double broadcast(double in,int rank_from=0);
#ifdef USE_MPI
  MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos);
  MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos);
  void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type);
#else
  uint64_t ceil_to_next_eight_multiple(uint64_t pos);
  uint64_t diff_with_next_eight_multiple(uint64_t pos);
#endif
  
  inline bool is_master_rank()
  {
    return rank==master_rank;
  }
  
  std::string MPI_get_processor_name();
  
  /// Component used to store a rank id
  struct MpiRank :
    BaseComp<MpiRank,int,0>
  {
    using Base=BaseComp<MpiRank,int,0>;
    using Base::Base;
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Send and receive a trivial type
  template <TriviallyCopyable T>
  T mpiSendrecv(const size_t& rankTo,
		const T& toSend,
		const size_t& rankFr)
  {
    T toRecv;
    
    MPI_Sendrecv(&toSend,sizeof(T),MPI_CHAR,rankTo,0,
		 &toRecv,sizeof(T),MPI_CHAR,rankFr,0,
		 MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    return toRecv;
  }
  
  /// Specialization of sendrecv for std::vector
  template <TriviallyCopyable T>
  std::vector<T> mpiSendrecv(const size_t& rankTo,
			     const std::vector<T>& toSend,
			     const size_t& rankFr)
  {
    const size_t nToSend=
      toSend.size();
    
    const size_t nToRecv=
      mpiSendrecv(rankTo,nToSend,rankFr);
    
    std::vector<T> toRecv(nToRecv);
    MPI_Sendrecv(&toSend[0],nToSend*sizeof(T),MPI_CHAR,rankTo,0,
		 &toRecv[0],nToRecv*sizeof(T),MPI_CHAR,rankFr,0,
		 MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    return toRecv;
  }
}

#endif
