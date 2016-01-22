#ifndef _IOS_HPP
#define _IOS_HPP

#include <string>
#include <stdio.h>
#include <stdint.h>
#include "new_types/complex.hpp"

#ifndef EXTERN_IOS
 #define EXTERN_IOS extern
#endif

namespace nissa
{
  EXTERN_IOS int verb_call;
  EXTERN_IOS int verbosity_lv;
  
  int count_substrings(const char *str,const char *sub);
  FILE* open_file(std::string path,const char *mode);
  FILE* open_text_file_for_output(std::string path);
  int cd(std::string path);
  int cp(std::string path_out,std::string path_in);
  int create_dir(std::string path);
  int master_fprintf(FILE *stream,const char *format,...);
  //int rm(const char *path);
  std::string combine(const char *format,...);
  void close_file(FILE *file);
  void fprintf_friendly_filesize(FILE *fout,uint64_t quant);
  void fprintf_friendly_units(FILE *fout,uint64_t quant,uint64_t orders,const char *units);
  void take_last_characters(char *out,const char *in,int size);
  int count_file_lines(std::string path);
  int get_file_size(std::string path);
  void print_contraction_to_file(FILE *fout,int op1,int op2,complex *contr,int twall,const char *tag,double norm);
  void print_contractions_to_file(FILE *fout,int ncontr,const int *op1,const int *op2,complex *contr,int twall,const char *tag,double norm);
}

#undef EXTERN_IOS

#endif
