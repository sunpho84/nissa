#pragma once

//to be moved elsewhere soon

double sqr(double a)
{return a*a;}

int min_int(int a,int b)
{if(a<b) return a;else return b;}
int max_int(int a,int b)
{if(a>b) return a;else return b;}

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
  
  if(rank==0) fprintf(fout,"%d %s%s",quant,units_prefix[iord],units);
}

void master_printf_box(const char *template,...)
{
  if(rank==0)
    {
      char temp_out[1024];
      va_list ap;
      va_start(ap,template);
      vsnprintf(temp_out,1024,template,ap);
      va_end(ap);
      int l=strlen(temp_out);

      printf("\n /");
      for(int i=0;i<l;i++) printf("-");
      printf("\\\n |%s|\n \\",temp_out);
      for(int i=0;i<l;i++) printf("-");
      printf("/\n");
    }
}

//create a dir
int create_dir(char *path)
{
  char command[1024];
  sprintf(command,"mkdir -p %s",path);  
  int res=(rank==0) ? system(command) : 0;
  MPI_Bcast(&res,1,MPI_INT,0,MPI_COMM_WORLD);
  
  return res;
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

//Open a text file for output
FILE* open_text_file_for_output(char *outfile)
{
  FILE *fout=NULL;
  if(rank==0) fout=fopen(outfile,"w");
  if(rank==0 && fout==NULL)
    {
      fprintf(stderr,"Couldn't open the file: %s",outfile);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  return fout;
}

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
  char *buf=nissa_malloc("buf",sel,char);
  
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
  for(int i=1;i<=n;i++) log_factn+=logf(i);
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
