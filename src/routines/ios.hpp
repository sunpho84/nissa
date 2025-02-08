#ifndef _IOS_HPP
#define _IOS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>
#include <fstream>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "base/debug.hpp"
#include "base/randomDevice.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "routines/mpi_routines.hpp"

#if THREADS_TYPE == OPENMP_THREADS
# include <omp.h>
#endif

#ifndef EXTERN_IOS
# define EXTERN_IOS extern
#endif

#define RED_HIGHLIGHT				\
  "\x1b[31m"

#define YELLOW_HIGHLIGHT			\
  "\x1b[33m"

#define GREEN_HIGHLIGHT				\
  "\x1b[32m"

#define DO_NOT_HIGHLIGHT			\
  "\x1b[0m"

#define MASTER_PRINTF(ARGS...)			\
  ::nissa::master_fprintf(stdout,ARGS)

//wrappers for verbosity_lv?
#if MAX_VERBOSITY_LV>=1
# define VERBOSITY_LV1_MASTER_PRINTF(...) MACRO_GUARD(if(VERBOSITY_LV1) MASTER_PRINTF(__VA_ARGS__);)
#else
# define VERBOSITY_LV1_MASTER_PRINTF(...) MACRO_GUARD()
#endif

#if MAX_VERBOSITY_LV>=2
# define VERBOSITY_LV2_MASTER_PRINTF(...) MACRO_GUARD(if(VERBOSITY_LV2) MASTER_PRINTF(__VA_ARGS__);)
#else
# define VERBOSITY_LV2_MASTER_PRINTF(...) MACRO_GUARD()
#endif

#if MAX_VERBOSITY_LV>=3
# define VERBOSITY_LV3_MASTER_PRINTF(...) MACRO_GUARD(if(VERBOSITY_LV3) MASTER_PRINTF(__VA_ARGS__);)
#else
# define VERBOSITY_LV3_MASTER_PRINTF(...) MACRO_GUARD()
#endif

namespace nissa
{
  extern int rank;
  
  EXTERN_IOS int prepend_time;
  EXTERN_IOS int verb_call;
  
  int count_substrings(const char *str,const char *sub);
  FILE* open_file(std::string path,const char *mode,int ext_rank=master_rank);
  FILE* open_text_file_for_output(std::string path);
  int cd(std::string path);
  int cp(std::string path_out,std::string path_in);
  int create_dir(std::string path);
  __attribute__((format (printf,2,3)))
  int master_fprintf(FILE *stream,const char *format,...);
  //int rm(const char *path);
  __attribute__((format (printf,1,2)))
  std::string combine(const char *format,...);
  void master_get_temp_file(FILE *&fout,std::string &prefix);
  void close_file(FILE *file);
  void fprintf_friendly_filesize(FILE *fout,uint64_t quant);
  void fprintf_friendly_units(FILE *fout,uint64_t quant,uint64_t orders,const char *units);
  void take_last_characters(char *out,const char *in,int size);
  int count_file_lines(std::string path);
  int get_file_size(std::string path);
  void print_contraction_to_file(FILE *fout,int gso,int gsi,complex *contr,int twall,const char *tag,double norm,int skip_header=false);
  void print_contractions_to_file(FILE *fout,int ncontr,const int *gso,const int *gsi,complex *contr,int twall,const char *tag,double norm,int skip_header=false);
  inline void print_contractions_to_file(FILE *fout,std::vector<idirac_pair_t> &list,complex *contr,int twall,const char *tag,double norm,int skip_header=false)
  {
    int ncontr=list.size();
    int gso[ncontr],gsi[ncontr];
    for(int i=0;i<ncontr;i++)
      {
	gso[i]=list[i].so;
	gsi[i]=list[i].si;
      }
    print_contractions_to_file(fout,ncontr,gso,gsi,contr,twall,tag,norm,skip_header);
  }
  
  //read from a file, opened only on master rank
  template <class T> T master_fscan(FILE *stream,const char *tag)
  {
    int thread_id=
#if THREADS_TYPE == OPENMP_THREADS
      omp_get_thread_num()
#else
      0
#endif
      ;
    
    //scan
    T out=0;
    if(is_master_rank() and thread_id==0)
      if(fscanf(stream,tag,&out)!=1)
	CRASH("Unable to read!");
    
    //broadcast
    MPI_Bcast(&out,sizeof(T),MPI_CHAR,master_rank,MPI_COMM_WORLD);
    
    return out;
  }
  
  //read using a path
  template <class T> T master_fscan(std::string path,const char *tag)
  {
    FILE *stream=open_file(path,"r");
    T out=master_fscan<T>(stream,tag);
    close_file(stream);
    return out;
  }
  
  //read an integer with either a path or a file
  template <class T> int master_fscan_int(T par)
  {
    return master_fscan<int>(par,"%d");
  }
  
  //read a double with either a path or a file
  template <class T> double master_fscan_double(T par)
  {
    return master_fscan<double>(par,"%lg");
  }
  
  //create and check lock files
  template <class T=uint64_t>
  class lock_file_t
  {
    //store whether is inited
    bool inited{false};
    
    //variable containing the lock word
    T tag;
    
    //path to lock
    std::string path;
    
    void assert_inited()
    {
      if(not inited) CRASH("Needs to be inited");
    }
    
  public:
    
    //create the tag
    void init()
    {
      MASTER_PRINTF("Initializing the tag for a %zu bytes lock-file\n",sizeof(T));
      get_system_random(tag);
      
      inited=true;
    }
    
    //try to open and write the tag
    bool try_lock(const std::string &ext_path)
    {
      assert_inited();
      
      //store the lock path
      path=ext_path;
      
      //create the lock on master
      int written=true;
      if(is_master_rank())
	if(std::ofstream(path)<<tag<<std::endl)
	  written=true;
	else
	  written=false;
      else
	written=true;
      
      MPI_Bcast(&written,1,MPI_INT,master_rank,MPI_COMM_WORLD);
      
      if(written) MASTER_PRINTF("Created lock file %s\n",path.c_str());
      else MASTER_PRINTF("Failed to create the lock file %s\n",path.c_str());
      
      return written;
    }
    
    //try to open and read the tag
    bool check_lock()
    {
      assert_inited();
      //temporary tag read
      T test_tag;
      memset(&test_tag,0,sizeof(T));
      
      //read on master
      if(is_master_rank()) std::ifstream(path)>>test_tag;
      
      //broadcast
      MPI_Bcast(&test_tag,sizeof(T),MPI_CHAR,master_rank,MPI_COMM_WORLD);
      
      //return the comparison
      return (test_tag==tag);
    }
  };
  
  //avoid c++ warning
  template <typename...Args>
  void safe_snprintf(char *buf,int n,const char *pattern,const Args&...args)
  {
    if(snprintf(buf,n,pattern,args...)<0) CRASH("witing to %d long array",n);
  }
}

#undef EXTERN_IOS

#endif
