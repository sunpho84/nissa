#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#include <stdio.h>
#include <stdint.h>

#ifndef EXTERN_VECTORS
 #define EXTERN_VECTORS extern
#endif

//vector tags name
#define BORDERS_ALLOCATED 1
#define BORDERS_VALID 2
#define EDGES_ALLOCATED 4
#define EDGES_VALID 8
#define BORDERS_COMMUNICATED_AT_LEAST_ONCE 16
#define DO_NOT_SET_FLAGS 1
#define SEND_BACKWARD_BORD 1
#define SEND_FORWARD_BORD 2

#define NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED 1

#define NISSA_VECT_STRING_LENGTH 20
#if defined(BGQ) && !defined(BGQ_EMU)
 #define NISSA_VECT_ALIGNMENT 32
#else
 #define NISSA_VECT_ALIGNMENT 16
#endif

#define nissa_malloc(a,b,c) (c*)internal_nissa_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define nissa_free(a) internal_nissa_free((char**)&(a),__FILE__,__LINE__)

#define CRASH_IF_NOT_ALIGNED(a,b) MACRO_GUARD(if((long long int)(void*)a%b!=0) crash("alignement problem");)
#define IF_VECT_NOT_INITIALIZED() if(main_arr!=((char*)&main_vect)+sizeof(nissa_vect))

namespace nissa
{
  //nissa vector
  struct nissa_vect
  {
    int64_t nel;
    int64_t size_per_el;
    
    char tag[NISSA_VECT_STRING_LENGTH];
    char type[NISSA_VECT_STRING_LENGTH];
    
    char file[NISSA_VECT_STRING_LENGTH];
    int line;
    
    nissa_vect *prev;
    nissa_vect *next;
    
    uint32_t flag;
    
    //padding to keep memory alignment
    char pad[(NISSA_VECT_ALIGNMENT-(2*sizeof(int64_t)+3*NISSA_VECT_STRING_LENGTH+sizeof(int)+2*sizeof(nissa_vect*)+sizeof(uint32_t))%NISSA_VECT_ALIGNMENT)%
	      NISSA_VECT_ALIGNMENT];
  };
  
  EXTERN_VECTORS int warn_if_not_disallocated;
  //vectors
  EXTERN_VECTORS int64_t max_required_memory;
  EXTERN_VECTORS int64_t required_memory;
  EXTERN_VECTORS void *main_arr;
  EXTERN_VECTORS nissa_vect main_vect;
  EXTERN_VECTORS nissa_vect *last_vect;
  EXTERN_VECTORS void *return_malloc_ptr;
  
  char *get_vect_name(void *v);
  int check_borders_allocated(void *data);
  int check_borders_communicated_at_least_once(void *data);
  int check_borders_valid(void *data);
  int check_edges_allocated(void *data);
  int check_edges_valid(void *data);
  int64_t compute_vect_memory_usage();
  int get_vect_flag(void *v,unsigned int flag);
  nissa_vect* get_vect(void *v);
  void *internal_nissa_malloc(const char *tag,int64_t nel,int64_t size_per_el,const char *type,const char *file,int line);
  void crash_if_borders_not_allocated(void *v);
  void crash_if_edges_not_allocated(void *v);
  void ignore_borders_communications_warning(void *data);
  void initialize_main_vect();
  void internal_nissa_free(char **arr,const char *file,int line);
  void vector_copy(void *a,const void *b);
  void vector_reset(void *a);
  void last_vect_content_printf();
  void vect_content_fprintf(FILE *fout,nissa_vect *vect);
  void vect_content_printf(nissa_vect *vect);
  void print_all_vect_content();
  void reorder_vector(char *vect,int *order,int nel,int sel);
  void set_borders_invalid(void *data);
  void set_borders_valid(void *data);
  void set_edges_invalid(void *data);
  void set_edges_valid(void *data);
  void set_vect_flag(void *v,unsigned int flag);
  void set_vect_flag_non_blocking(void *v,unsigned int flag);
  void unset_vect_flag(void *v,unsigned int flag);
  void unset_vect_flag_non_blocking(void *v,unsigned int flag);
  void vect_content_fprintf(FILE *f,void *vec);
  void vect_content_printf(void *vec);
}

#undef EXTERN_VECTORS

#endif

