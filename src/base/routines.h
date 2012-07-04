#ifndef _ROUTINES_H
#define _ROUTINES_H
FILE* open_file(const char *outfile,const char *mode);
FILE* open_text_file_for_output(const char *outfile);
MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos);
MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos);
double glb_reduce_double(double in_loc);
double lfact(double n);
double max_double(double a,double b);
double min_double(double a,double b);
double sqr(double a);
double take_time();
int create_dir(char *path);
int factorize(int *list,int N);
int log2N(int N);
int master_fprintf(FILE *stream,const char *format,...);
int max_int(int a,int b);
int min_int(int a,int b);
void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type);
void fprintf_friendly_filesize(FILE *fout,int quant);
void fprintf_friendly_units(FILE *fout,int quant,int orders,const char *units);
void master_printf_box(const char *templ,...);
void reorder_vector(char *vect,int *order,int nel,int sel);
void swap_doubles(double *d1,double *d2);
void take_last_characters(char *out,const char *in,int size);
#endif
