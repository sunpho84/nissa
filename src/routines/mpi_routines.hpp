#ifndef _MPI_ROUTINES_HPP
#define _MPI_ROUTINES_HPP

#include <mpi.h>
#include <algorithm>

#include <routines/math_routines.hpp>
#include <new_types/coords.hpp>
#include <new_types/float_128.hpp>
#include <routines/rank.hpp>

#ifndef EXTERN_MPI
 #define EXTERN_MPI extern
 #define INIT_MPI_TO(...)
#else
 #define INIT_MPI_TO(ARGS...) ARGS
#endif

namespace nissa
{
  struct rat_approx_t;
  
  DECLARE_COMPONENT(Rank,int,DYNAMIC);
  
  /// Coordinates of a rank
  using RankCoords=Coords<Rank>;
  
  //ranks
  EXTERN_MPI RankCoords fix_nranks;
  CUDA_MANAGED EXTERN_MPI RankCoords rank_coord;
  EXTERN_MPI RankCoords rank_neigh[2],rank_neighdw,rank_neighup;
  EXTERN_MPI RankCoords plan_rank,line_rank,line_coord_rank;
  CUDA_MANAGED EXTERN_MPI RankCoords nrank_dir;
  //basic mpi types
  EXTERN_MPI MPI_Datatype MPI_FLOAT_128;
  EXTERN_MPI MPI_Datatype MPI_COMPLEX_128;
  EXTERN_MPI MPI_Datatype MPI_SU3;
  EXTERN_MPI MPI_Datatype MPI_QUAD_SU3;
  EXTERN_MPI MPI_Datatype MPI_AS2T_SU3;
  EXTERN_MPI MPI_Datatype MPI_COLOR;
  EXTERN_MPI MPI_Datatype MPI_SPIN;
  EXTERN_MPI MPI_Datatype MPI_SPINSPIN;
  EXTERN_MPI MPI_Datatype MPI_SPINCOLOR;
  EXTERN_MPI MPI_Datatype MPI_SPINCOLOR_128;
  EXTERN_MPI MPI_Datatype MPI_REDSPINCOLOR;
  //float 128 summ
  EXTERN_MPI MPI_Op MPI_FLOAT_128_SUM;
  EXTERN_MPI MPI_Op MPI_COMPLEX_128_SUM;
  
  EXTERN_MPI MPI_Datatype MPI_LX_SU3_EDGES_SEND[NDIM*(NDIM-1)/2],MPI_LX_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_LX_AS2T_SU3_EDGES_SEND[NDIM*(NDIM-1)/2],MPI_LX_AS2T_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_LX_QUAD_SU3_EDGES_SEND[NDIM*(NDIM-1)/2],MPI_LX_QUAD_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_EO_QUAD_SU3_EDGES_SEND[96],MPI_EO_QUAD_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  EXTERN_MPI MPI_Datatype MPI_EO_AS2T_SU3_EDGES_SEND[96],MPI_EO_AS2T_SU3_EDGES_RECE[NDIM*(NDIM-1)/2];
  
  //volume, plan and line communicator
  EXTERN_MPI MPI_Comm cart_comm;
  EXTERN_MPI MPI_Comm plan_comm[NDIM];
  EXTERN_MPI MPI_Comm line_comm[NDIM];
  
#define DEFINE_MPI_DATATYPE_OF(T,MPI_T)		\
  /*! MPI Datatype corresponding to T */	\
  inline MPI_Datatype _MPI_Datatype_of(T*)	\
  {						\
    return MPI_T;				\
  }
  
  DEFINE_MPI_DATATYPE_OF(int64_t,MPI_LONG)
  DEFINE_MPI_DATATYPE_OF(double,MPI_DOUBLE)
  DEFINE_MPI_DATATYPE_OF(complex,MPI_DOUBLE_COMPLEX)
  DEFINE_MPI_DATATYPE_OF(float_128,MPI_FLOAT_128)
  DEFINE_MPI_DATATYPE_OF(complex_128,MPI_COMPLEX_128)
  
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
  
  DEFINE_MPI_OP_DISPATCHER(float_128,MPI_FLOAT_128_SUM);
  DEFINE_MPI_OP_DISPATCHER(complex_128,MPI_COMPLEX_128_SUM);
  
  /// Gets the sum operation for the type T
  template <typename T>
  MPI_Op MPI_Op_sum_for_type()
  {
   return _MPI_Op_dispatcher<T>::sum();
  }
  
#undef DEFINE_MPI_OP_DISPATCHER
  
  size_t MPI_Get_count_size_t(MPI_Status &status);
  
  //broadcast a coord
  template <typename I>
  void coords_broadcast(Coords<I>& c)
  {
    FOR_ALL_DIRECTIONS(mu)
      {
	auto& ref=c(mu)();
	MPI_Bcast(&ref,NDIM,MPI_Datatype_of<std::decay_t<decltype(ref)>>(),master_rank,MPI_COMM_WORLD);
      }
  }
  
  void get_MPI_nranks();
  void get_MPI_rank();
  void init_MPI_thread(int narg,char **arg);
  void define_MPI_types();
  void create_MPI_cartesian_grid();
  void ranks_abort(int err);
  void ranks_barrier();
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
  
  std::string MPI_get_processor_name();
}

#endif
