//to be moved elsewhere soon

#include <mpi.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string>

#include "debug.h"
#include "global_variables.h"
#include "random.h"
#include "routines.h"
#include "vectors.h"
#include "../new_types/float128.h"
#include "../new_types/new_types_definitions.h"

double sqr(double a)
{return a*a;}

int min_int(int a,int b)
{if(a<b) return a;else return b;}
int max_int(int a,int b)
{if(a>b) return a;else return b;}

//take the different with following multiple of eight
MPI_Offset diff_with_next_eight_multiple(MPI_Offset pos)
{
  MPI_Offset diff=pos%8;
  if(diff!=0) diff=8-diff;

  return diff;
}

//perform the metropolis test on passed value
int metro_test(double arg)
{
  double tresh=(arg<=0) ? 1 : exp(-arg);
  double toss=rnd_get_unif(&glb_rnd_gen,0,1);
  
  return toss<tresh;
}

//ceil to next multiple of eight
MPI_Offset ceil_to_next_eight_multiple(MPI_Offset pos)
{return pos+diff_with_next_eight_multiple(pos);}

void MPI_FLOAT_128_SUM_routine(void *in,void *out,int *len,MPI_Datatype *type)
{for(int i=0;i<(*len);i++) float_128_summassign(((float_128*)out)[i],((float_128*)in)[i]);}

double min_double(double a,double b)
{if(a<b) return a;else return b;}
double max_double(double a,double b)
{if(a>b) return a;else return b;}

//print a number in some kind of familiar units
void fprintf_friendly_units(FILE *fout,int quant,int orders,const char *units)
{
  const char units_prefix[6][2]={"","K","M","G","T","P"};
  int iord=0;
  double temp=quant;
  
  while(temp>orders)
    {
      temp/=orders;
      iord++;
    }
  
  quant=(int)(temp+0.5);
  
  master_fprintf(fout,"%d %s%s",quant,units_prefix[iord],units);
}

void master_printf_box(const char *templ,...)
{
  char temp_out[1024];
  va_list ap;
  va_start(ap,templ);
  vsnprintf(temp_out,1024,templ,ap);
  va_end(ap);
  int l=strlen(temp_out);
  
  master_printf("\n /");
  for(int i=0;i<l;i++) master_printf("-");
  printf("\\\n |%s|\n \\",temp_out);
  for(int i=0;i<l;i++) master_printf("-");
  master_printf("/\n");
}

//create a dir
int create_dir(char *path)
{
  int res=(rank==0) ? mkdir(path,480) : 0;
  MPI_Bcast(&res,1,MPI_INT,0,MPI_COMM_WORLD);
  if(res!=0)
    master_printf("Warning, failed to create dir %s, returned %d. Check that you have permissions and that parent dir exists.\n",path,res);
  
  return res;
}

//broadcast an int
int master_broadcast(int in)
{
  MPI_Bcast(&in,1,MPI_INT,0,MPI_COMM_WORLD);
  return in;
}

//copy a file
int cp(char *path_out,char *path_in)
{
  int rc;
  if(rank==0)
    {
      char command[1024];
      sprintf(command,"cp %s %s",path_in,path_out);
      rc=system(command);
    }
  
  return master_broadcast(rc);
}

//remove a file
int rm(const char *path)
{
  int rc;
  if(rank==0)
    {
      char command[1024];
      sprintf(command,"rm %s",path);
      rc=system(command);
    }
  
  return master_broadcast(rc);
}

//pass to the folder
int cd(const char *path)
{
  int rc;
  if(rank==0)
    {
      char command[1024];
      sprintf(command,"cd %s",path);
      rc=system(command);
    }
  
  return master_broadcast(rc);
}

void fprintf_friendly_filesize(FILE *fout,int quant)
{fprintf_friendly_units(fout,quant,1024,"Bytes");}

//swap two doubles
void swap_doubles(double *d1,double *d2)
{
  double tmp=(*d1);
  (*d1)=(*d2);
  (*d2)=tmp;
}

double take_time()
{
  //MPI_Barrier(MPI_COMM_WORLD);
  return MPI_Wtime();
}

//Open a file checking it
FILE* open_file(const char *outfile,const char *mode)
{
  FILE *fout=NULL;
  
  if(rank==0)
    {
      fout=fopen(outfile,mode);
      if(fout==NULL) crash("Couldn't open the file: %s for mode: %s",outfile,mode);
    }
  
  return fout;
}

//Open a text file for output
FILE* open_text_file_for_output(const char *outfile)
{return open_file(outfile,"w");}

void take_last_characters(char *out,const char *in,int size)
{
  int len=strlen(in)+1;
  int copy_len=(len<=size)?len:size;
  const char *str_init=(len<=size)?in:in+len-size;
  memcpy(out,str_init,copy_len);
}

//reorder a vector according to the specified order (the order is destroyed)
void reorder_vector(char *vect,int *order,int nel,int sel)
{
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

int master_fprintf(FILE *stream,const char *format,...)
{
  int ret=0;
  
  if(rank==0)
    {
      va_list ap;
      va_start(ap,format);
      ret=vfprintf(stream,format,ap);
      va_end(ap);
    }
  
  return ret;
}

//return the log of the factorial of n
double lfact(double n)
{
  double log_factn=0;
  for(int i=1;i<=n;i++) log_factn+=log(i);
  return log_factn;
}

//return the log2 of N
int log2N(int N)
{
  int log2N=0;
  
  do log2N++;
  while ((2<<log2N)<N);
  
  return log2N;
}

//factorize a number
int factorize(int *list,int N)
{
  int nfatt=0;
  int fatt=2;

  while(N>1)
    {
      int div=N/fatt;
      int res=N-div*fatt;
      if(res!=0) fatt++;
      else 
	{
	  N=div;
	  list[nfatt]=fatt;
	  nfatt++;
	}
    }

  return nfatt;
}

//reduce a complex
void glb_reduce_complex(complex out_glb,complex in_loc)
{MPI_Allreduce(in_loc,out_glb,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}

//reduce a double
double glb_reduce_double(double in_loc)
{
  double out_glb;
  
  MPI_Allreduce(&in_loc,&out_glb,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return out_glb;
}

//reduce a double vector
void glb_reduce_double_vect(double *out_glb,double *in_loc,int nel)
{MPI_Allreduce(in_loc,out_glb,nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}

//reduce a complex vector
void glb_reduce_complex_vect(complex *out_glb,complex *in_loc,int nel)
{MPI_Allreduce((double*)in_loc,(double*)out_glb,2*nel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}

//reduce an int
int glb_reduce_int(int in_loc)
{
  int out_glb;
  
  MPI_Allreduce(&in_loc,&out_glb,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  return out_glb;
}

void set_gauge_action_type(theory_pars_type &theory_pars,char *type)
{
  if(strcmp(type,"Wilson")==0) theory_pars.gac_type=Wilson_action;
  else
    if(strcmp(type,"tlSym")==0) theory_pars.gac_type=tlSym_action;
    else crash("unknown gauge action: %s",type);
}

std::string combine(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  return std::string(buffer);
}
