#ifndef _IOS_HPP
#define _IOS_HPP

#include <string>
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <mpi.h>

#include "new_types/complex.hpp"
#include "operations/contract/mesons_2pts.hpp"

#ifdef USE_THREADS
 #include <omp.h>
#endif

#ifndef EXTERN_IOS
 #define EXTERN_IOS extern
#endif

#define master_printf(...) nissa::master_fprintf(stdout,__VA_ARGS__)

//add verbosity macro
#if MAX_VERBOSITY_LV>=1
 #define VERBOSITY_LV1 (nissa::verbosity_lv>=1)
#else
 #define VERBOSITY_LV1 0
#endif
#if MAX_VERBOSITY_LV>=2
 #define VERBOSITY_LV2 (nissa::verbosity_lv>=2)
#else
 #define VERBOSITY_LV2 0
#endif
#if MAX_VERBOSITY_LV>=3
 #define VERBOSITY_LV3 (nissa::verbosity_lv>=3)
#else
 #define VERBOSITY_LV3 0
#endif

#define NISSA_DEFAULT_VERBOSITY_LV 1

//wrappers for verbosity_lv?
#define verbosity_lv1_master_printf(...) MACRO_GUARD(if(VERBOSITY_LV1) master_printf(__VA_ARGS__);)
#define verbosity_lv2_master_printf(...) MACRO_GUARD(if(VERBOSITY_LV2) master_printf(__VA_ARGS__);)
#define verbosity_lv3_master_printf(...) MACRO_GUARD(if(VERBOSITY_LV3) master_printf(__VA_ARGS__);)

namespace nissa
{
  extern int rank;
  
  EXTERN_IOS int verb_call;
  EXTERN_IOS int verbosity_lv;
  
  int count_substrings(const char *str,const char *sub);
  FILE* open_file(std::string path,const char *mode,int ext_rank=0);
  FILE* open_text_file_for_output(std::string path);
  int cd(std::string path);
  int cp(std::string path_out,std::string path_in);
  int create_dir(std::string path);
  int master_fprintf(FILE *stream,const char *format,...);
  //int rm(const char *path);
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
#ifdef USE_THREADS
      omp_get_thread_num()
#else
      0
#endif
      ;
    
    //scan
    T out=0;
    if(rank==0 and thread_id==0)
      if(fscanf(stream,tag,&out)!=1)
	crash("Unable to read!");
    
    //broadcast
    MPI_Bcast(&out,sizeof(T),MPI_CHAR,0,MPI_COMM_WORLD);
    
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
  template <class T> int master_fscan_double(T par)
  {
    return master_fscan<int>(par,"%lg");
  }
}

#undef EXTERN_IOS

#endif
