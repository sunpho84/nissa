#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "new_types/float128.h"
#include "new_types/new_types_definitions.h"
#include "routines/ios.h"
#include "routines/math_routines.h"
#ifdef USE_THREADS
 #include "routines/thread.h"
#endif

#include "debug.h"
#include "global_variables.h"
#include "thread_macros.h"

//#define DEBUG

//return the pointer to the nissa vect
nissa_vect* get_nissa_vec(void *v)
{return (nissa_vect*)v-1;}

//return the name of the vector
char *get_vec_name(void *v)
{return get_nissa_vec(v)->tag;}

//print the content of an nissa vect
void nissa_vect_content_fprintf(FILE *fout,nissa_vect *vect)
{
  if(rank==0)
    {
      fprintf(fout,"\"%s\" ",vect->tag);
      fprintf(fout,"of %d elements of type \"%s\" (%d bytes) ",vect->nel,vect->type,vect->nel*vect->size_per_el);
      fprintf(fout,"allocated in file %s line %d\n",vect->file,vect->line);
    }
}
//wrappers
void nissa_vect_content_printf(nissa_vect *vect)
{nissa_vect_content_fprintf(stdout,vect);}
void vec_content_fprintf(FILE *f,void *vec)
{nissa_vect_content_fprintf(f,(nissa_vect*)vec-1);}
void vec_content_printf(void *vec)
{vec_content_fprintf(stdout,vec);}
void last_nissa_vect_content_printf()
{
  master_printf("Last nissa vect content: ");
  nissa_vect_content_printf(last_nissa_vect);
}

//set a flag: we need to be sure that all the threads are consistent
void set_vec_flag_non_blocking(void *v,unsigned int flag)
{get_nissa_vec(v)->flag|=flag;}
void set_vec_flag(void *v,unsigned int flag)
{
#ifdef DEBUG
  GET_THREAD_ID();
  printf("set_vec_flag for vect %s allocated in file %s line %d, rank %d thread_id: %d, thread_pool_locked: %d\n",
	 get_vec_name(v),get_nissa_vec(v)->file,get_nissa_vec(v)->line,rank,thread_id,thread_pool_locked);
  print_backtrace_list();
#endif
  //update atomically
  if((get_nissa_vec(v)->flag & flag)!=flag)
    {
      THREAD_BARRIER();
      get_nissa_vec(v)->flag|=flag;
    }
  THREAD_BARRIER();
}

//unset a flag (idem)
void unset_vec_flag_non_blocking(void *v,unsigned int flag)
{get_nissa_vec(v)->flag &= ~flag;}
void unset_vec_flag(void *v,unsigned int flag)
{
#ifdef DEBUG
  printf("unset_vec_flag for vect %s allocated in file %s line %d, rank %d thread_id: %d, thread_pool_locked: %d\n",
	 get_vec_name(v),get_nissa_vec(v)->file,get_nissa_vec(v)->line,rank,thread_id,thread_pool_locked);
#endif
  //update atomically
  if(((~get_nissa_vec(v)->flag)&flag)!=flag)
    {
      THREAD_BARRIER();
      get_nissa_vec(v)->flag&=~flag;
    }
  THREAD_BARRIER();
}

//get a flag
int get_vec_flag(void *v,unsigned int flag)
{return get_nissa_vec(v)->flag & flag;}

//check if edges are valid
int check_edges_valid(void *data)
{return get_vec_flag(data,EDGES_VALID);}

//check if edges have been allocated
int check_edges_allocated(void *data)
{return get_vec_flag(data,EDGES_ALLOCATED);}

//check if borders are valid
int check_borders_valid(void *data)
{return get_vec_flag(data,BORDERS_VALID);}

//check if borders have been allocated
int check_borders_allocated(void *data)
{return get_vec_flag(data,BORDERS_ALLOCATED);}

//check if borders have been communicated at least once
int check_borders_communicated_at_least_once(void *data)
{return get_vec_flag(data,BORDERS_COMMUNICATED_AT_LEAST_ONCE);}

//ignore the warning on borders allocated but never used
void ignore_borders_communications_warning(void *data)
{set_vec_flag(data,BORDERS_COMMUNICATED_AT_LEAST_ONCE);}

//set borders ad valid
void set_borders_valid(void *data)
{set_vec_flag(data,BORDERS_VALID|BORDERS_COMMUNICATED_AT_LEAST_ONCE);}

//set edges ad valid
void set_edges_valid(void *data)
{set_vec_flag(data,EDGES_VALID);}

//set borders as invalid
void set_borders_invalid(void *data)
{unset_vec_flag(data,BORDERS_VALID|EDGES_VALID);}

//set edges as invalid
void set_edges_invalid(void *data)
{unset_vec_flag(data,EDGES_VALID);}

//print all nissa vect
void print_all_nissa_vect_content()
{
  nissa_vect *curr=&(main_nissa_vect);
  do
    {  
      nissa_vect_content_printf(curr);
      curr=curr->next;
    }
  while(curr!=NULL);
}

//check if the borders are allocated
void crash_if_borders_not_allocated(void *v)
{
  if(!check_borders_allocated(v))
    if(rank==0)
      {
	fprintf(stderr,"borders not allocated in ");
	vec_content_fprintf(stderr,v);
	crash("see error");
      }
}

//check if the edges are allocated
void crash_if_edges_not_allocated(void *v)
{
  if(!check_edges_allocated(v))
    if(rank==0)
      {
	fprintf(stderr,"edges not allocated in ");
	vec_content_fprintf(stderr,v);
	crash("see error");
      }
}

//print all nissa vect
int compute_nissa_vect_memory_usage()
{
  int tot=0;
  nissa_vect *curr=&(main_nissa_vect);
  do
    {  
      tot+=curr->nel*curr->size_per_el;
      curr=curr->next;
    }
  while(curr!=NULL);
  
  return tot;
}

//initialize the first vector
void initialize_main_nissa_vect()
{
  if_nissa_vect_not_initialized()
    {
      nissa_required_memory=0;
      nissa_max_required_memory=0;
      last_nissa_vect=&main_nissa_vect;
      sprintf(main_nissa_vect.tag,"base");
      sprintf(main_nissa_vect.type,"(null)");
      main_nissa_vect.prev=main_nissa_vect.next=NULL;
      main_nissa_vect.nel=0;
      main_nissa_vect.size_per_el=0;
      memcpy(main_nissa_vect.file,__FILE__+std::max(0,(int)strlen(__FILE__)-12),12);
      main_nissa_vect.line=__LINE__;
      main_nissa_arr=(char*)last_nissa_vect+sizeof(nissa_vect);

      master_printf("Vector memory manager started\n");
    }
}

//allocate an nissa vector 
void *internal_nissa_malloc(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line)
{
  GET_THREAD_ID();
  if(IS_MASTER_THREAD)
    {
      if_nissa_vect_not_initialized() initialize_main_nissa_vect();
      
      int size=nel*size_per_el;
      //try to allocate the new vector
      nissa_vect *nv=(nissa_vect*)malloc(size+sizeof(nissa_vect));
      if(nv==NULL)
	crash("could not allocate vector named \"%s\" of %d elements of type %s (total size: %d bytes) "
	      "request on line %d of file %s",tag,nel,type,size,line,file);
      
      //fill the vector with information supplied
      nv->line=line;
      nv->nel=nel;
      nv->size_per_el=size_per_el;
      nv->flag=0;
      take_last_characters(nv->file,file,nissa_vect_string_length);
      take_last_characters(nv->tag,tag,nissa_vect_string_length);
      take_last_characters(nv->type,type,nissa_vect_string_length);
      
      //append the new vector to the list
      nv->next=NULL;
      nv->prev=last_nissa_vect;
      
      last_nissa_vect->next=nv;
      last_nissa_vect=nv;
      
      if(nissa_verbosity>=3)
	{
	  master_printf("Allocated vector ");
	  nissa_vect_content_printf(last_nissa_vect);
	}
      
      //define returned pointer and check for its alignement
      return_nissa_malloc_ptr=(void*)(last_nissa_vect+1);
      int offset=((long long int)(return_nissa_malloc_ptr))%nissa_vect_alignment;
      if(offset!=0)
	crash("memory alignment problem, vector %s has %d offset",tag,offset);
      
      //if borders or edges are allocated, set appropriate flag
      if(nel==(loc_vol+bord_vol) || nel==(loc_volh+bord_volh))
	set_vec_flag_non_blocking(return_nissa_malloc_ptr,BORDERS_ALLOCATED);
      if(nel==(loc_vol+bord_vol+edge_vol) || nel==(loc_volh+bord_volh+edge_volh))
	set_vec_flag_non_blocking(return_nissa_malloc_ptr,BORDERS_ALLOCATED|EDGES_ALLOCATED);
      
      //Update the amount of required memory
      nissa_required_memory+=size;
      nissa_max_required_memory=std::max(nissa_max_required_memory,nissa_required_memory);

      cache_flush();
    }
  
  //sync so we are sure that master thread allocated
  THREAD_BARRIER();
  void *res=return_nissa_malloc_ptr;
  
  //resync so all threads return the same pointer
  THREAD_BARRIER();
  
  return res;
}

//copy one vector into another
void vector_copy(void *a,void *b)
{
  GET_THREAD_ID();
  
  //sync so we are sure that all threads are here
  THREAD_BARRIER();
  
  if(IS_MASTER_THREAD)
    {
      
      //control that a!=b
      if(a==b) crash("while copying, vector1 and vector2 are the same");
      
      //get nissa-vector pointers
      nissa_vect *nissa_a=(nissa_vect*)((char*)a-sizeof(nissa_vect));
      nissa_vect *nissa_b=(nissa_vect*)((char*)b-sizeof(nissa_vect));
      
      //get nentries
      int nel_a=nissa_a->nel;
      int nel_b=nissa_b->nel;
      
      //get nentries per site
      int size_per_el_a=nissa_a->size_per_el;
      int size_per_el_b=nissa_b->size_per_el;
      
      //check size agreement
      if(nel_a!=nel_b) crash("while copying, vector %s contains %d and vector %s contains %d",nissa_a->tag,nel_a,nissa_b->tag,nel_b);
      
      //check type agreement
      if(size_per_el_a!=size_per_el_b)
	crash("while copying, vector %s contains %d bytes per el and vector %s contains %d",nissa_a->tag,size_per_el_a,nissa_b->tag,size_per_el_b);
      
      //perform the copy
      memcpy(a,b,size_per_el_a*nel_a);
      
      //copy the flag
      nissa_a->flag=nissa_b->flag;
    }

  //sync so we are sure that all threads are here
  THREAD_BARRIER();      
}

//reset a vector
void vector_reset(void *a)
{
  //sync so all thread are not using the vector
  THREAD_BARRIER();
  
  GET_THREAD_ID();
  if(IS_MASTER_THREAD)
    {
      nissa_vect *nissa_a=(nissa_vect*)((char*)a-sizeof(nissa_vect));
      int nel_a=nissa_a->nel;
      int size_per_el_a=nissa_a->size_per_el;
      memset(a,0,size_per_el_a*nel_a);
    }
  
  //sync so all thread see that have been reset
  THREAD_BARRIER();
}

//release a vector
void internal_nissa_free(char **arr,const char *file,int line)
{
  //sync so all thread are not using the vector
  THREAD_BARRIER();

  GET_THREAD_ID();
  if(IS_MASTER_THREAD)
    {
      if(arr!=NULL)
	{
	  nissa_vect *vect=(nissa_vect*)((char*)(*arr)-sizeof(nissa_vect));
	  nissa_vect *prev=vect->prev;
	  nissa_vect *next=vect->next;
	  
	  if(nissa_verbosity>=3)
	    {
	      master_printf("At line %d of file %s freeing vector ",line,file);
	      nissa_vect_content_printf(vect);
	    }
	  
	  if(nissa_warn_if_not_communicated)
	    if(nissa_nranks>1 && check_borders_allocated(*arr) && !check_borders_communicated_at_least_once(*arr))
	      master_printf("Warning, you allocated borders for vector: %s on line %d of file %s, but never communicated them!\n",vect->tag,vect->line,vect->file);
	  
	  //detach from previous
	  prev->next=next;
	  
	  //if not last element
	  if(next!=NULL) next->prev=prev;
	  else last_nissa_vect=prev;
	  
	  //update the nissa required memory
	  nissa_required_memory-=(vect->size_per_el*vect->nel);
	  
	  //really free
	  free(vect);
	}
      else crash("Error, trying to delocate a NULL vector on line: %d of file: %s\n",line,file);
      
      *arr=NULL;
      cache_flush();
    }
  
  //sync so all thread see that have deallocated
  THREAD_BARRIER();
}

//reorder a vector according to the specified order (the order is destroyed)
void reorder_vector(char *vect,int *order,int nel,int sel)
{
  GET_THREAD_ID();
  if(!IS_MASTER_THREAD) crash("not threaded yet");
  char *buf=(char*)nissa_malloc("buf",sel,char);
  
  for(int sour=0;sour<nel;sour++)
    while(sour!=order[sour])
      {
        int dest=order[sour];
        
        memcpy(buf,vect+sour*sel,sel);
        memcpy(vect+sour*sel,vect+dest*sel,sel);
        memcpy(vect+dest*sel,buf,sel);
        
        order[sour]=order[dest];
        order[dest]=dest;
      }
  
  set_borders_invalid(vect);
  
  nissa_free(buf);
}
