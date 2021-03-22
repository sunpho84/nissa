#ifndef _MPI_ROUTINES_HPP
#define _MPI_ROUTINES_HPP

#include <mpi.h>
#include <algorithm>

#include "geometry/geometry_lx.hpp"
#include "math_routines.hpp"
#include "new_types/float_128.hpp"

#ifndef EXTERN_MPI
 #define EXTERN_MPI extern
 #define INIT_MPI_TO(...)
#else
 #define INIT_MPI_TO(ARGS...) ARGS
#endif

namespace nissa
{
  struct rat_approx_t;
  
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
  
  EXTERN_MPI int master_rank INIT_MPI_TO(=0);
  
#define DEFINE_MPI_DATATYPE_OF(T,MPI_T)		\
  /*! MPI Datatype corresponding to T */	\
  inline MPI_Datatype _MPI_Datatype_of(T*)	\
  {						\
    return MPI_T;				\
  }
  
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
    static MPI_Op sum()				\
    {						\
      return OP;				\
    }						\
  }
  
  DEFINE_MPI_OP_DISPATCHER(float_128,MPI_FLOAT_128_SUM);
  DEFINE_MPI_OP_DISPATCHER(complex_128,MPI_COMPLEX_128_SUM);
  
#undef DEFINE_MPI_OP_DISPATCHER
  
  size_t MPI_Get_count_size_t(MPI_Status &status);
  void coords_broadcast(coords c);
  void get_MPI_nranks();
  void get_MPI_rank();
  void init_MPI_thread(int narg,char **arg);
  void define_MPI_types();
  void create_MPI_cartesian_grid();
  void ranks_abort(int err);
  void ranks_barrier();
  float glb_reduce_single(float in_loc);
  double glb_reduce_double(double in_loc,double (*thread_op)(const double&,const double&)=summ<double>,MPI_Op mpi_op=MPI_SUM);
  inline double glb_max_double(double in_loc)
  {return glb_reduce_double(in_loc,nissa_max,MPI_MAX);}
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
  
  inline bool is_master_rank()
  {
    return rank==master_rank;
  }
  
  std::string MPI_get_processor_name();
  
  /////////////////////////////////////////////////////////////////
  
  //reduce a double
  inline float MPI_reduce_single(float in_loc)
  {
    float out_glb;
    
    MPI_Allreduce(&in_loc,&out_glb,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce an int
  inline int MPI_reduce_int(int in_loc)
  {
    int out_glb;
    
    MPI_Allreduce(&in_loc,&out_glb,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce a vector
  template <typename T>
  inline void MPI_reduce_vect(T* out_glb,T* in_loc=nullptr,const int n=1,MPI_Op mpi_op=MPI_SUM)
  {
    if(in_loc==nullptr)
      in_loc=(T*)MPI_IN_PLACE;
    
    MPI_Allreduce(in_loc,out_glb,n,MPI_Datatype_of<T>(),mpi_op,MPI_COMM_WORLD);
  }
  
  //reduce
  template <typename T>
  inline void MPI_reduce(T* out_glb,T* in_loc,MPI_Op mpi_op=_MPI_Op_dispatcher<T>::sum())
  {
    MPI_reduce_vect(out_glb,in_loc,1,mpi_op);
  }
  
  //reduce
  template <typename T>
  inline void MPI_reduce(T* out_glb,MPI_Op mpi_op=_MPI_Op_dispatcher<T>::sum())
  {
    MPI_reduce_vect(out_glb,nullptr,1,mpi_op);
  }
  
  void glb_nodes_reduce_double_vect(double *out_glb,double *in_loc,int nel);
  inline void glb_nodes_reduce_double_vect(double *vect,int nel)
  {glb_nodes_reduce_double_vect(vect,(double*)MPI_IN_PLACE,nel);}
  inline void glb_nodes_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel)
  {glb_nodes_reduce_double_vect(out_glb[0],in_loc[0],2*nel);}
  inline void glb_nodes_reduce_complex_vect(complex *vect,int nel)
  {
    glb_nodes_reduce_double_vect(vect[0],2*nel);
  }
  void glb_reduce_float_128(float_128 out_glb,float_128 in_loc);
  
  void glb_reduce_complex_128(complex_128 out_glb,complex_128 in_loc);
}

#endif
