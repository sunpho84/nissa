#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cwctype>
#include <stdio.h>
#include <stdint.h>

#include "debug.hpp"

#ifndef EXTERN_VECTORS
# define EXTERN_VECTORS extern
#endif

//vector tags name
#define DO_NOT_SET_FLAGS 1

#define NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED 1

#define NISSA_VECT_STRING_LENGTH 20
#define NISSA_VECT_ALIGNMENT 16

#define nissa_malloc(a,b,c) (c*)internal_nissa_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define nissa_free(a) do{static_assert(std::is_pointer_v<std::remove_reference_t<decltype((a))>>,"What are you deallocating?");internal_nissa_free((char**)&(a),__FILE__,__LINE__);}while(0)

#define CRASH_IF_NOT_ALIGNED(a,b) MACRO_GUARD(if((long long int)(void*)a%b!=0) CRASH("alignement problem");)
#define IF_MAIN_VECT_NOT_INITIALIZED() if(main_arr!=((char*)&main_vect)+sizeof(nissa_vect))

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
    
    void assert_is_nissa_vect() const
    {
      const auto check=
	[](const char* v)
	{
	  if(not iswalpha(v[0]))
	     return false;
	  
	  int i=1;
	  
	  bool f=false;
	  while((i<NISSA_VECT_STRING_LENGTH) and not f)
	    if(v[i++]=='\0')
	      f=true;
	  
	  return f;
	};
      
      if(not (check(tag) and check(type)))
	CRASH("not a nissa_vector");
    }
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
  int64_t compute_vect_memory_usage();
  int get_vect_flag(void *v,unsigned int flag);
  nissa_vect* get_vect(void *v);
  void *internal_nissa_malloc(const char *tag,int64_t nel,int64_t size_per_el,const char *type,const char *file,int line);
  void initialize_main_vect();
  void internal_nissa_free(char **arr,const char *file,int line);
  void vector_copy(void *a,const void *b);
  void vector_reset(void *a);
  void last_vect_content_printf();
  void vect_content_fprintf(FILE *fout,nissa_vect *vect);
  void vect_content_printf(nissa_vect *vect);
  void print_all_vect_content();
  void set_vect_flag(void *v,unsigned int flag);
  void set_vect_flag_non_blocking(void *v,unsigned int flag);
  void unset_vect_flag(void *v,unsigned int flag);
  void unset_vect_flag_non_blocking(void *v,unsigned int flag);
  void vect_content_fprintf(FILE *f,void *vec);
  void vect_content_printf(void *vec);
}

#undef EXTERN_VECTORS

#endif

