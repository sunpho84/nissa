#pragma once

//to be moved elsewhere soon

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
  MPI_Barrier(MPI_COMM_WORLD);
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
