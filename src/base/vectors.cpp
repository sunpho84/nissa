#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "debug.hpp"

#define EXTERN_VECTORS
 #include "vectors.hpp"

//#define DEBUG

namespace nissa
{
  //return the pointer to the nissa vect
  nissa_vect* get_vect(void *v)
  {
    return (nissa_vect*)v-1;
  }
  
  //return the name of the vector
  char *get_vect_name(void *v)
  {
    return get_vect(v)->tag;
  }
  
  //print the content of an nissa vect
  void vect_content_fprintf(FILE *fout,nissa_vect *vect)
  {
    if(is_master_rank())
      {
	fprintf(fout,"\"%s\" ",vect->tag);
	fprintf(fout,"of %ld elements of type \"%s\" (%ld bytes) ",vect->nel,vect->type,vect->nel*vect->size_per_el);
	fprintf(fout,"allocated in file %s line %d\n",vect->file,vect->line);
      }
  }
  //wrappers
  void vect_content_printf(nissa_vect *vect)
  {vect_content_fprintf(stdout,vect);}
  void vect_content_fprintf(FILE *f,void *vec)
  {vect_content_fprintf(f,(nissa_vect*)vec-1);}
  void vect_content_printf(void *vec)
  {vect_content_fprintf(stdout,vec);}
  void last_vect_content_printf()
  {
    master_printf("Last vect content: ");
    vect_content_printf(last_vect);
  }
  
  //set a flag: we need to be sure that all the threads are consistent
  void set_vect_flag_non_blocking(void *v,unsigned int flag)
  {get_vect(v)->flag|=flag;}
  void set_vect_flag(void *v,unsigned int flag)
  {
#ifdef DEBUG
    printf("set_vect_flag for vect %s allocated in file %s line %d, rank %d thread_id: %d, thread_pool_locked: %d\n",
	   get_vect_name(v),get_vect(v)->file,get_vect(v)->line,rank,thread_id,thread_pool_locked);
    print_backtrace_list();
#endif
    //update atomically
    if((get_vect(v)->flag & flag)!=flag)
      {
	//THREAD_BARRIER();
	get_vect(v)->flag|=flag;
      }
    //THREAD_BARRIER();
  }
  
  //unset a flag (idem)
  void unset_vect_flag_non_blocking(void *v,unsigned int flag)
  {get_vect(v)->flag &= ~flag;}
  void unset_vect_flag(void *v,unsigned int flag)
  {
#ifdef DEBUG
    printf("unset_vect_flag for vect %s allocated in file %s line %d, rank %d thread_id: %d, thread_pool_locked: %d\n",
	   get_vect_name(v),get_vect(v)->file,get_vect(v)->line,rank,thread_id,thread_pool_locked);
#endif
    //update atomically
    if(((~get_vect(v)->flag)&flag)!=flag)
      {
	//THREAD_BARRIER();
	get_vect(v)->flag&=~flag;
      }
    //THREAD_BARRIER();
  }
  
  //get a flag
  int get_vect_flag(void *v,unsigned int flag)
  {return get_vect(v)->flag & flag;}
  
  //check if edges are valid
  int check_edges_valid(void *data)
  {return get_vect_flag(data,EDGES_VALID);}
  
  //check if edges have been allocated
  int check_edges_allocated(void *data,int min_size)
  {
    return ((get_vect(data)->nel*get_vect(data)->size_per_el)>=min_size) or
      get_vect_flag(data,EDGES_ALLOCATED);
  }
  
  //check if borders are valid
  int check_borders_valid(void *data)
  {return get_vect_flag(data,BORDERS_VALID);}
  
  //check if borders have been allocated
  int check_borders_allocated(void *data,int min_size)
  {
    return ((get_vect(data)->nel*get_vect(data)->size_per_el)>=min_size) or
      get_vect_flag(data,BORDERS_ALLOCATED);
  }
  
  //check if borders have been communicated at least once
  int check_borders_communicated_at_least_once(void *data)
  {return get_vect_flag(data,BORDERS_COMMUNICATED_AT_LEAST_ONCE);}
  
  //ignore the warning on borders allocated but never used
  void ignore_borders_communications_warning(void *data)
  {set_vect_flag(data,BORDERS_COMMUNICATED_AT_LEAST_ONCE);}
  
  //set borders ad valid
  void set_borders_valid(void *data)
  {set_vect_flag(data,BORDERS_VALID|BORDERS_COMMUNICATED_AT_LEAST_ONCE);}
  
  //set edges ad valid
  void set_edges_valid(void *data)
  {set_vect_flag(data,EDGES_VALID);}
  
  //set borders as invalid
  void set_borders_invalid(void *data)
  {unset_vect_flag(data,BORDERS_VALID|EDGES_VALID);}
  
  //set edges as invalid
  void set_edges_invalid(void *data)
  {unset_vect_flag(data,EDGES_VALID);}
  
  //print all nissa vect
  void print_all_vect_content()
  {
    nissa_vect *curr=&(main_vect);
    do
      {
	vect_content_printf(curr);
	curr=curr->next;
      }
    while(curr!=NULL);
  }
  
  //check if the borders are allocated
  void crash_if_borders_not_allocated(void *v,int min_size)
  {
    if(!check_borders_allocated(v,min_size))
      if(is_master_rank())
	{
	  fprintf(stderr,"borders not allocated in ");
	  vect_content_fprintf(stderr,v);
	  crash("see error");
	}
  }
  
  //check if the edges are allocated
  void crash_if_edges_not_allocated(void *v,int min_size)
  {
    if(!check_edges_allocated(v,min_size))
      if(is_master_rank())
	{
	  fprintf(stderr,"edges not allocated in ");
	  vect_content_fprintf(stderr,v);
	  crash("see error");
	}
  }
  
  //print all nissa vect
  int64_t compute_vect_memory_usage()
  {
    int64_t tot=0;
    nissa_vect *curr=&(main_vect);
    do
      {
	tot+=curr->nel*curr->size_per_el;
	curr=curr->next;
      }
    while(curr!=NULL);
    
    return tot;
  }
  
  //initialize the first vector
  void initialize_main_vect()
  {
    IF_MAIN_VECT_NOT_INITIALIZED()
      {
	required_memory=0;
	max_required_memory=0;
	last_vect=&main_vect;
	sprintf(main_vect.tag,"base");
	sprintf(main_vect.type,"(null)");
	main_vect.prev=main_vect.next=NULL;
	main_vect.nel=0;
	main_vect.size_per_el=0;
	const char file_name[]=__FILE__;
	memcpy(main_vect.file,file_name+std::max(0,(int)strlen(file_name)-12),12);
	main_vect.line=__LINE__;
	main_arr=(char*)last_vect+sizeof(nissa_vect);
	
	master_printf("Vector memory manager started\n");
      }
  }
  
  //allocate an nissa vector
  void *internal_nissa_malloc(const char *tag,int64_t nel,int64_t size_per_el,const char *type,const char *file,int line)
  {
    if(IS_MASTER_THREAD)
      {
	IF_MAIN_VECT_NOT_INITIALIZED() initialize_main_vect();
	
	if(VERBOSITY_LV3)
	  {
	    master_printf("Allocating vector ");
	    vect_content_printf(last_vect);
	  }
	
	int64_t size=nel*size_per_el;
	//try to allocate the new vector
	nissa_vect *nv;
	int64_t tot_size=size+sizeof(nissa_vect);
#define ALLOCATING_ERROR \
	"could not allocate vector named \"%s\" of %d elements of type %s (total size: %d bytes) "\
		"request on line %d of file %s",tag,nel,type,size,line,file
#if THREADS_TYPE==CUDA_THREADS
	decript_cuda_error(cudaMallocManaged(&nv,tot_size),ALLOCATING_ERROR);
#else
	nv=(nissa_vect*)malloc(tot_size);
	if(nv==NULL)
	  crash(ALLOCATING_ERROR);
#endif
#undef ALLOCATING_ERROR
	
	//fill the vector with information supplied
	nv->line=line;
	nv->nel=nel;
	nv->size_per_el=size_per_el;
	nv->flag=0;
	take_last_characters(nv->file,file,NISSA_VECT_STRING_LENGTH);
	take_last_characters(nv->tag,tag,NISSA_VECT_STRING_LENGTH);
	take_last_characters(nv->type,type,NISSA_VECT_STRING_LENGTH);
	
	//append the new vector to the list
	nv->next=NULL;
	nv->prev=last_vect;
	
	last_vect->next=nv;
	last_vect=nv;
	
	//define returned pointer and check for its alignement
	return_malloc_ptr=(void*)(last_vect+1);
	int64_t offset=((int64_t)(return_malloc_ptr))%NISSA_VECT_ALIGNMENT;
	if(offset!=0)
	  crash("memory alignment problem, vector %s has %d offset",tag,offset);
	
	//if borders or edges are allocated, set appropriate flag
	if(nel==(locVol+bord_vol) || nel==(locVolh+bord_volh))
	  set_vect_flag_non_blocking(return_malloc_ptr,BORDERS_ALLOCATED);
	if(nel==(locVol+bord_vol+edge_vol) || nel==(locVolh+bord_volh+edge_volh))
	  set_vect_flag_non_blocking(return_malloc_ptr,BORDERS_ALLOCATED|EDGES_ALLOCATED);
	
	//Update the amount of required memory
	required_memory+=size;
	max_required_memory=std::max(max_required_memory,required_memory);
	
	//cache_flush();
      }
    
    //sync so we are sure that master thread allocated
    //THREAD_BARRIER();
    void *res=return_malloc_ptr;
    
    //resync so all threads return the same pointer
    //THREAD_BARRIER();
    
    return res;
  }
  
  //copy one vector into another
  void vector_copy(void *a,const void *b)
  {
    
    //sync so we are sure that all threads are here
    //THREAD_BARRIER();
    
    //control that a!=b
    if(a!=b)
      {
	//get vector pointers
	nissa_vect *nissa_a=(nissa_vect*)((char*)a-sizeof(nissa_vect));
	nissa_vect *nissa_b=(nissa_vect*)((char*)b-sizeof(nissa_vect));
	
	//get nentries
	int64_t nel_a=nissa_a->nel;
	int64_t nel_b=nissa_b->nel;
	
	//get nentries per site
	int64_t size_per_el_a=nissa_a->size_per_el;
	int64_t size_per_el_b=nissa_b->size_per_el;
	
	//check size agreement
	if(nel_a!=nel_b) crash("while copying, vector %s allocated at line %d of file %s contains %d and vector %s allocated at line %d of file %s contains %d",
			       nissa_a->tag,nissa_a->line,nissa_a->file,nel_a,nissa_b->tag,nissa_b->line,nissa_b->file,nel_b);
	
	//check type agreement
	if(size_per_el_a!=size_per_el_b)
	  crash("while copying, vector %s contains %d bytes per el and vector %s contains %d",
		nissa_a->tag,size_per_el_a,nissa_b->tag,size_per_el_b);
	
	//perform the copy
	NISSA_PARALLEL_LOOP(i,0,nel_a)
	  memcpy((char*)a+i*size_per_el_a,(char*)b+i*size_per_el_a,size_per_el_a);
	NISSA_PARALLEL_LOOP_END;
	
	//copy the flag
	nissa_a->flag=nissa_b->flag;
	
	//sync so we are sure that all threads are here
	//THREAD_BARRIER();
      }
  }
  
  //reset a vector
  void vector_reset(void *a)
  {
    if(IS_MASTER_THREAD)
      {
	const nissa_vect* nissa_a=(nissa_vect*)((char*)a-sizeof(nissa_vect));
	nissa_a->assert_is_nissa_vect();
	
	const int64_t nel_a=nissa_a->nel;
	const int64_t size_per_el_a=nissa_a->size_per_el;
	memset(a,0,size_per_el_a*nel_a);
      }
  }
  
  //release a vector
  void internal_nissa_free(char **arr,const char *file,int line)
  {
    if(arr!=NULL)
      {
	nissa_vect *vect=get_vect(*arr);
	nissa_vect *prev=vect->prev;
	nissa_vect *next=vect->next;
	
	if(VERBOSITY_LV3)
	  {
	    master_printf("At line %d of file %s freeing vector ",line,file);
	    vect_content_printf(vect);
	  }
	
	if(warn_if_not_communicated)
	  if(nranks>1 && check_borders_allocated(*arr,0) && !check_borders_communicated_at_least_once(*arr))
	    master_printf("Warning, you allocated borders for vector: %s on line %d of file %s, but never communicated them!\n",vect->tag,vect->line,vect->file);
	
	//detach from previous
	prev->next=next;
	
	//if not last element
	if(next!=NULL) next->prev=prev;
	else last_vect=prev;
	
	//update the required memory
	required_memory-=(vect->size_per_el*vect->nel);
	
	//really free
#if THREADS_TYPE == CUDA_THREADS
	decript_cuda_error(cudaFree(vect),"freeing the memory for vector");
#else
	free(vect);
#endif
      }
    else crash("Error, trying to delocate a NULL vector on line: %d of file: %s\n",line,file);
    
    *arr=NULL;
  }
  
  //reorder a vector according to the specified order
  void reorder_vector(char* vect,int* order,int nel,int sel)
  {
    char *buf=nissa_malloc("buf",(int64_t)sel*nel,char);
    
    //copy in the buffer
    NISSA_PARALLEL_LOOP(sour,0,nel)
      memcpy(buf+order[sour]*sel,vect+sour*sel,sel);
    NISSA_PARALLEL_LOOP_END;
    
    //THREAD_BARRIER();
    
    NISSA_PARALLEL_LOOP(sour,0,nel)
      memcpy(vect+sour*sel,buf+sour*sel,sel);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(vect);
    
    nissa_free(buf);
  }
  
}
