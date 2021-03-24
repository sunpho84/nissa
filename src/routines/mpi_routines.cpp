#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#define EXTERN_MPI
 #include "mpi_routines.hpp"

#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"
#include "new_types/rat_approx.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //return  the count covnerted to size_t
  size_t MPI_Get_count_size_t(MPI_Status &status)
  {
    int nbytes;
    decript_MPI_error(MPI_Get_count(&status,MPI_BYTE,&nbytes),"while counting bytes");
    if(nbytes<0) crash("negative count: %d",nbytes);
    
    return (size_t)nbytes;
  }
  
  //take the different with following multiple of eight
  MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos)
  {
    MPI_Offset diff=pos%8;
    if(diff!=0) diff=8-diff;
    
    return diff;
  }
  
  //summ two float_128
  void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type)
  {for(int i=0;i<(*len);i++) float_128_summassign(((float_128*)out)[i],((float_128*)in)[i]);}
  
  //summ two complex_128
  void MPI_COMPLEX_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type)
  {for(int i=0;i<(*len);i++) complex_128_summassign(((complex_128*)out)[i],((complex_128*)in)[i]);}
  
  //init mpi
  void init_MPI_thread(int narg,char **arg)
  {
#ifdef USE_MPI
    
 #if THREADS_TYPE != NO_THREADS
    int provided;
    MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
 #else
    MPI_Init(&narg,&arg);
 #endif
#endif
  }
  
  //get nranks
  void get_MPI_nranks()
  {
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nranks);
#else
    nranks=1;
#endif
  }
  
  //get rank
  void get_MPI_rank()
  {
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
    rank=0;
#endif
  }
  
  //define the cartesian grid
  void create_MPI_cartesian_grid()
  {
#ifdef USE_MPI
    coords periods;
    for(int mu=0;mu<NDIM;mu++) periods[mu]=1;
    MPI_Cart_create(MPI_COMM_WORLD,NDIM,nrank_dir,periods,1,&cart_comm);
    //takes rank and ccord of local rank
    MPI_Comm_rank(cart_comm,&cart_rank);
    MPI_Cart_coords(cart_comm,cart_rank,NDIM,rank_coord);
    
    //create communicator along plan
    for(int mu=0;mu<NDIM;mu++)
      {
	coords split_plan;
	coords proj_rank_coord;
	for(int nu=0;nu<NDIM;nu++)
	  {
	    split_plan[nu]=(nu==mu) ? 0 : 1;
	    proj_rank_coord[nu]=(nu==mu) ? 0 : rank_coord[nu];
	  }
	MPI_Cart_sub(cart_comm,split_plan,&(plan_comm[mu]));
	MPI_Comm_rank(plan_comm[mu],&(plan_rank[mu]));
	if(plan_rank[mu]!=rank_of_coord(proj_rank_coord))
	  crash("Plan communicator has messed up coord: %d and rank %d (implement reorder!)",
		rank_of_coord(proj_rank_coord),plan_rank[mu]);
      }
    
    //create communicator along line
    for(int mu=0;mu<NDIM;mu++)
      {
	//split the communicator
	coords split_line;
	memset(split_line,0,sizeof(coords));
	split_line[mu]=1;
	MPI_Cart_sub(cart_comm,split_line,&(line_comm[mu]));
	
	//get rank id
	MPI_Comm_rank(line_comm[mu],&(line_rank[mu]));
	
	//get rank coord along line comm
	MPI_Cart_coords(line_comm[mu],line_rank[mu],1,&(line_coord_rank[mu]));
	
	//check communicator
	if(line_rank[mu]!=rank_coord[mu] || line_rank[mu]!=line_coord_rank[mu])
	  crash("Line communicator has messed up coord and rank (implement reorder!)");
      }
#else
    cart_rank=plan_rank=line_rank=0;
    for(int mu=0;mu<NDIM;mu++) rank_coord[mu]=planline_coord[mu]=0;
#endif
  }
  
  //barrier
  void ranks_barrier()
  {
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  //abort
  void ranks_abort(int err)
  {
#ifdef USE_MPI
    printf("thread %d on rank %d aborting\n",THREAD_ID,rank);
    MPI_Abort(MPI_COMM_WORLD,0);
#else
    exit(0);
#endif
  }
  //define all types
  void define_MPI_types()
  {
#ifdef USE_MPI
    //128 bit float
    MPI_Type_contiguous(2,MPI_DOUBLE,&MPI_FLOAT_128);
    MPI_Type_commit(&MPI_FLOAT_128);
    
    //128 bit complex
    MPI_Type_contiguous(2,MPI_FLOAT_128,&MPI_COMPLEX_128);
    MPI_Type_commit(&MPI_COMPLEX_128);
    
    //define the gauge link
    MPI_Type_contiguous(18,MPI_DOUBLE,&MPI_SU3);
    MPI_Type_commit(&MPI_SU3);
    
    //four (in 4d) links starting from a single point
    MPI_Type_contiguous(NDIM,MPI_SU3,&MPI_QUAD_SU3);
    MPI_Type_commit(&MPI_QUAD_SU3);
    
    //six (in 4d) links starting from a single point
    MPI_Type_contiguous(sizeof(as2t_su3)/sizeof(su3),MPI_SU3,&MPI_AS2T_SU3);
    MPI_Type_commit(&MPI_AS2T_SU3);
    
    //a color (6 doubles)
    MPI_Type_contiguous(6,MPI_DOUBLE,&MPI_COLOR);
    MPI_Type_commit(&MPI_COLOR);
    
    //a spin (8 doubles)
    MPI_Type_contiguous(8,MPI_DOUBLE,&MPI_SPIN);
    MPI_Type_commit(&MPI_SPIN);
    
    //a spinspin (32 doubles)
    MPI_Type_contiguous(32,MPI_DOUBLE,&MPI_SPINSPIN);
    MPI_Type_commit(&MPI_SPINSPIN);
    
    //a spincolor (24 doubles)
    MPI_Type_contiguous(24,MPI_DOUBLE,&MPI_SPINCOLOR);
    MPI_Type_commit(&MPI_SPINCOLOR);
    
    //a spincolor_128 (24 float_128)
    MPI_Type_contiguous(24,MPI_FLOAT_128,&MPI_SPINCOLOR_128);
    MPI_Type_commit(&MPI_SPINCOLOR_128);
    
    //a reduced spincolor (12 doubles)
    MPI_Type_contiguous(12,MPI_DOUBLE,&MPI_REDSPINCOLOR);
    MPI_Type_commit(&MPI_REDSPINCOLOR);
    
    //summ for 128 bit float
    MPI_Op_create((MPI_User_function*)MPI_FLOAT_128_SUM_routine,1,&MPI_FLOAT_128_SUM);
    
    //summ for 128 bit complex
    MPI_Op_create((MPI_User_function*)MPI_COMPLEX_128_SUM_routine,1,&MPI_COMPLEX_128_SUM);
#endif
  }
  
  //broadcast a coord
  void coords_broadcast(coords c)
  {MPI_Bcast(c,NDIM,MPI_INT,master_rank,MPI_COMM_WORLD);}
  
  //ceil to next multiple of eight
  MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos)
  {return pos+diff_with_next_eight_multiple(pos);}
  
  //internal version
  template <class T> T broadcast_internal(T in,int rank_from,MPI_Datatype type)
  {
    T out;
    if(IS_MASTER_THREAD) MPI_Bcast(&in,1,type,rank_from,MPI_COMM_WORLD);
    THREAD_BROADCAST(out,in);
    return out;
  }
  
  //broadcast an int
  int broadcast(int in,int rank_from)
  {return broadcast_internal(in,rank_from,MPI_INT);}
  
  //broadcast an int
  double broadcast(double in,int rank_from)
  {return broadcast_internal(in,rank_from,MPI_DOUBLE);}
  
  //broadcast a whole rational approximation
  void broadcast(rat_approx_t *rat,int rank_from)
  {
    
    //first destroy on non-sending
    THREAD_BARRIER(); //need to barrier to avoid race condition when later "rat-degree" is update through mpi_bcast
    
    //get degree
    if(IS_MASTER_THREAD)
      {
	int degree=0;
	
	if(rank_from==rank) degree=rat->degree();
	MPI_Bcast(&degree,1,MPI_INT,rank_from,MPI_COMM_WORLD);
	
	//allocate if not generated here
	if(rank_from!=rank)	rat->resize(degree);
	
	//and now broadcast the remaining part
	MPI_Bcast(rat->name,20,MPI_CHAR,rank_from,MPI_COMM_WORLD);
	MPI_Bcast(&rat->minimum,1,MPI_DOUBLE,rank_from,MPI_COMM_WORLD);
	MPI_Bcast(&rat->maximum,1,MPI_DOUBLE,rank_from,MPI_COMM_WORLD);
    	MPI_Bcast(&rat->maxerr,1,MPI_DOUBLE,rank_from,MPI_COMM_WORLD);
    	MPI_Bcast(&rat->num,1,MPI_INT,rank_from,MPI_COMM_WORLD);
    	MPI_Bcast(&rat->den,1,MPI_INT,rank_from,MPI_COMM_WORLD);
    	MPI_Bcast(&rat->cons,1,MPI_DOUBLE,rank_from,MPI_COMM_WORLD);
	MPI_Bcast(rat->poles.data(),rat->degree(),MPI_DOUBLE,rank_from,MPI_COMM_WORLD);
	MPI_Bcast(rat->weights.data(),rat->degree(),MPI_DOUBLE,rank_from,MPI_COMM_WORLD);
      }
    THREAD_BARRIER();
  }
  
  //Return the name of the processor
  std::string MPI_get_processor_name()
  {
    int resultlen=MPI_MAX_PROCESSOR_NAME;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name,&resultlen);
    
    return name;
  }
}
