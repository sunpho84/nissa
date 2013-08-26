#ifndef _INPUT_H
#define _INPUT_H

int dir_exists(char *path);
int file_exists(const char *path);
int read_var_catcherr(char *out,const char *par,int size_of);
void close_input();
void expect_str(const char *exp_str);
void file_touch(const char *path);
void open_input(char *input_path);
void read_double(double *out);
void read_int(int *out);
void read_list_of_chars(const char *tag,int *nentries,char ***list,int nchar_per_entry);
void read_list_of_double_pairs(const char *tag,int *nentries,double **list1,double **list2);
void read_list_of_doubles(const char *tag,int *nentries,double **list);
void read_list_of_ints(const char *tag,int *nentries,int **list);
void read_list_of_var(const char *tag,int *nentries,char **list,int size_of_el,const char *par);
void read_list_of_var_pairs(const char *tag,int *nentries,char **list1,char **list2,int size_of_el,const char *par);
void read_nissa_config_file();
void read_str(char *str,int length);
void read_str_double(const char *exp_str,double *in);
void read_str_momentum_t(const char *exp_str,momentum_t in);
void read_str_int(const char *exp_str,int *in);
void read_str_str(const char *exp_str,char *in,int length);
void read_var(char *out,const char *par,int size_of);

#endif
