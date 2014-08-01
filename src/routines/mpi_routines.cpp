#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#include "ios.hpp"

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"
#include "new_types/new_types_definitions.hpp"

#ifdef USE_THREADS
  #include "routines/thread.hpp"
#endif

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
  
  //init mpi
  void init_MPI_thread(int narg,char **arg)
  {
#ifdef USE_MPI

 #ifdef USE_THREADS
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
    int periods[4]={1,1,1,1};
    MPI_Cart_create(MPI_COMM_WORLD,4,nrank_dir,periods,1,&cart_comm);
    //takes rank and ccord of local rank
    MPI_Comm_rank(cart_comm,&cart_rank);
    MPI_Cart_coords(cart_comm,cart_rank,4,rank_coord);
    
    //create communicator along plan
    for(int mu=0;mu<4;mu++)
      {
	int split_plan[4];
	coords proj_rank_coord;
	for(int nu=0;nu<4;nu++)
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
    for(int mu=0;mu<4;mu++)
      {
	//split the communicator
	int split_line[4];
	memset(split_line,0,4*sizeof(int));
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
    for(int mu=0;mu<4;mu++) rank_coord[mu]=planline_coord[mu]=0;
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
    GET_THREAD_ID();
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
    
    //define the gauge link
    MPI_Type_contiguous(18,MPI_DOUBLE,&MPI_SU3);
    MPI_Type_commit(&MPI_SU3);
    
    //four links starting from a single point
    MPI_Type_contiguous(4,MPI_SU3,&MPI_QUAD_SU3);
    MPI_Type_commit(&MPI_QUAD_SU3);
    
    //six links starting from a single point
    MPI_Type_contiguous(6,MPI_SU3,&MPI_AS2T_SU3);
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
#endif
  }
  
  //broadcast a coord
  void coords_broadcast(coords c)
  {MPI_Bcast(c,4,MPI_INT,0,MPI_COMM_WORLD);}
  
  //ceil to next multiple of eight
  MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos)
  {return pos+diff_with_next_eight_multiple(pos);}
  
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
    
#ifdef USE_THREADS
    if(!thread_pool_locked)
      {
	GET_THREAD_ID();
	
	//copy loc in the buf and sync all the threads
	glb_double_reduction_buf[thread_id]=in_loc;
	THREAD_BARRIER();
	
	//within master thread summ all the pieces and between MPI
	if(IS_MASTER_THREAD)
	  {
	    for(unsigned int ith=1;ith<nthreads;ith++) in_loc+=glb_double_reduction_buf[ith];
	    MPI_Allreduce(&in_loc,&(glb_double_reduction_buf[0]),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    cache_flush();
	  }
	
	//read glb val
	THREAD_ATOMIC_EXEC(out_glb=glb_double_reduction_buf[0];);
      }
    else
#endif
      MPI_Allreduce(&in_loc,&out_glb,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce a double
  float glb_reduce_single(float in_loc)
  {
    float out_glb;
    
#ifdef USE_THREADS
    if(!thread_pool_locked)
      {
	GET_THREAD_ID();
	
	//copy loc in the buf and sync all the threads
	glb_single_reduction_buf[thread_id]=in_loc;
	THREAD_BARRIER();
	
	//within master thread summ all the pieces and between MPI
	if(IS_MASTER_THREAD)
	  {
	    for(unsigned int ith=1;ith<nthreads;ith++) in_loc+=glb_single_reduction_buf[ith];
	    MPI_Allreduce(&in_loc,&(glb_single_reduction_buf[0]),1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	    cache_flush();
	  }
	
	//read glb val
	THREAD_ATOMIC_EXEC(out_glb=glb_single_reduction_buf[0];);
      }
    else
#endif
      MPI_Allreduce(&in_loc,&out_glb,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //max of all double
  double glb_max_double(double in_loc)
  {
    double out_glb;
    
#ifdef USE_THREADS
    if(!thread_pool_locked)
      {
	GET_THREAD_ID();
	
	//copy loc in the buf and sync all the threads
	glb_double_reduction_buf[thread_id]=in_loc;
	THREAD_BARRIER();
	
	//within master thread summ all the pieces and between MPI
	if(IS_MASTER_THREAD)
	  {
	    for(unsigned int ith=1;ith<nthreads;ith++) in_loc=std::max(in_loc,glb_double_reduction_buf[ith]);
	    MPI_Allreduce(&in_loc,&(glb_double_reduction_buf[0]),1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	    cache_flush();
	  }
	
	//read glb val
	THREAD_ATOMIC_EXEC(out_glb=glb_double_reduction_buf[0];);
      }
    else
#endif
      MPI_Allreduce(&in_loc,&out_glb,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    
    return out_glb;
  }
  
  //reduce an int
  void glb_reduce_int(int *out_glb,int in_loc)
  {
#ifdef USE_THREADS
    if(!thread_pool_locked) crash("not threaded yet");
    else
#endif
      MPI_Allreduce(&in_loc,out_glb,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  
  //reduce a complex
  void glb_reduce_complex(complex out_glb,complex in_loc)
  {for(int ri=0;ri<2;ri++) out_glb[ri]=glb_reduce_double(in_loc[ri]);}
  
  //reduce a float_128
  void glb_reduce_float_128(float_128 out_glb,float_128 in_loc)
  {
#ifdef USE_THREADS
    if(!thread_pool_locked)
      {
	GET_THREAD_ID();
	
	//copy loc in the buf and sync all the threads
	float_128_copy(glb_quadruple_reduction_buf[thread_id],in_loc);
	if(VERBOSITY_LV3 && rank==0) printf(" entering 128 reduction, in_loc[%d]: %+016.16lg\n",
					    THREAD_ID,in_loc[0]+in_loc[1]);
	THREAD_BARRIER();
	
	//within master thread summ all the pieces and between MPI
	if(IS_MASTER_THREAD)
	  {
	    for(unsigned int ith=1;ith<nthreads;ith++) float_128_summassign(in_loc,glb_quadruple_reduction_buf[ith]);
	    MPI_Allreduce(in_loc,glb_quadruple_reduction_buf[0],1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
	    if(VERBOSITY_LV3 && rank==0) printf("glb tot: %+016.16lg\n",glb_quadruple_reduction_buf[0][0]+
						glb_quadruple_reduction_buf[0][1]);
	    cache_flush();
	  }
	
	//read glb val
	THREAD_ATOMIC_EXEC(float_128_copy(out_glb,glb_quadruple_reduction_buf[0]));
      }
    else
#endif
      MPI_Allreduce(in_loc,out_glb,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
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
}
