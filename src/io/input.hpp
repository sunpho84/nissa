#ifndef _INPUT_HPP
#define _INPUT_HPP

#ifndef EXTERN_INPUT
 #define EXTERN_INPUT extern
 #define INIT_EXTERN_TO(A)
#else
 #define INIT_EXTERN_TO(A) =A
#endif

#include <stdio.h>
#include <string>
#include <vector>

#include "geometry/geometry_lx.hpp"

namespace nissa
{
  //global input file handle
  EXTERN_INPUT std::vector<FILE*> input_global_stack;
  EXTERN_INPUT FILE* input_global INIT_EXTERN_TO(NULL);
  
  int file_lock(std::string path);
  int file_unlock(int f);
  int dir_exists(std::string path);
  int fileExists(const std::string& path);
  int read_var_catcherr(char *out,const char *par,int size_of);
  void close_input();
  void expect_str(const char *exp_str);
  void file_touch(std::string path);
  void open_input(std::string input_path);
  void read_double(double *out);
  void read_int(int *out);
  void read_int64_t(int64_t *out);
  void read_list_of_chars(const char *tag,int *nentries,char ***list,int nchar_per_entry);
  void read_list_of_double_pairs(const char *tag,int *nentries,double **list1,double **list2);
  void read_list_of_doubles(const char *tag,int *nentries,double **list);
  void read_list_of_double_triples(const char *tag,int *nentries,double **list1,double **list2,double **list3);
  void read_list_of_ints(const char *tag,int *nentries,int **list);
  void read_list_of_var(const char *tag,int *nentries,char **list,int size_of_el,const char *par);
  void read_list_of_var_pairs(const char *tag,int *nentries,char **list1,char **list2,int size_of_el,const char *par);
  void read_list_of_var_triples(const char *tag,int *nentries,char **list1,char **list2,char **list3,int size_of_el,const char *par);
  void read_nissa_config_file();
  void read_str(char *str,int length);
  void read_str_double(const char *exp_str,double *in);
  void read_str_momentum_t(const char *exp_str,Momentum& in);
  void read_str_int(const char *exp_str,int *in);
  void read_str_str(const char *exp_str,char *in,int length);
  void read_var(char *out,const char *par,int size_of);
}

#endif
