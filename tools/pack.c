#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

typedef char name[10];
int T;

void crash(const char *templ,...)
{
  va_list ap;
  va_start(ap,templ);
  
  char mess[1024];
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"%s\n",mess);  
  
  exit(1);
}

FILE *open_file(const char *path,const char *mode)
{
  FILE *out=fopen(path,mode);
  if(out==NULL) crash("error opening file %s in mode %s",path,mode);
  
  return out;
}

FILE *open_pipe(const char *path,const char *mode)
{
  FILE *out=popen(path,mode);
  if(out==NULL) crash("error opening pipe with command %s in mode %s",path,mode);
  
  return out;
}

void copy_lines(FILE *fout,FILE *fin,int n)
{
  for(int i=0;i<n;i++)
    {
      char line[1024];
      fgets(line,1024,fin);
      
      fprintf(fout,"%s",line);
      //printf("%s",line);
    }
}

//count the number of combo
int count_combo(const char *path_in)
{
  char command[1024];
  sprintf(command,"grep -c m1 %s",path_in);
  FILE *pop=open_pipe(command,"r");
  
  int nc,nr=fscanf(pop,"%d",&nc);
  if(nr!=1) crash("error counting combo, returned %d from grep, %d from scanf",nr,nc);
  
  pclose(pop);
  
  return nc;
}

//fully parse the file
void parse_file(const char *path_out,const char *path_in,int ncorr,name *corr_name,int *npieces)
{
  int ncombo=count_combo(path_in);
  printf("File contains %d combo\n",ncombo);
  
  FILE *fout=open_file(path_out,"w");
  FILE *fin=open_file(path_in,"r");
  
  for(int icombo=0;icombo<ncombo;icombo++)
    {
      //header of the combo
      copy_lines(fout,fin,2);
      
      for(int icorr=0;icorr<ncorr;icorr++)
	{
	  double corr[T*2];
	  for(int i=0;i<T*2;i++) corr[i]=0;
	  
	  //print out header
	  fprintf(fout," # %s\n",corr_name[icorr]);
	  
	  for(int ipiece=0;ipiece<npieces[icorr];ipiece++)
	    {
	      //discard the name of the corr
	      char line[1024];
	      fgets(line,1024,fin);
	      //printf("discarded corr name: %s",line);
	      //read the corr piece
	      for(int i=0;i<T*2;i++)
		{
		  double buf;
		  int nr=fscanf(fin,"%lg",&buf);
		  if(nr!=1) crash("error reading, returned %d, read %lg",nr,buf);
		  
		  //printf("%lg\n",buf);
		  corr[i]+=buf/npieces[icorr];
		}
	      
	      //read upon end of line
	      fgets(line,1024,fin);
	      //printf("%s",line);
	      
	      //discard following empty line
	      fgets(line,1024,fin);
	      //printf("%s",line);
	    }
	  
	  //print corr
	  for(int t=0;t<T;t++) fprintf(fout,"%+016.16g\t%+016.16g\n",corr[t*2],corr[t*2+1]);
	  fprintf(fout,"\n");
	}
    }
  
  fclose(fin);
  fclose(fout);
}

int main(int narg,char **arg)
{
  if(narg<2) crash("error, use %s T",arg[0]);
  sscanf(arg[1],"%d",&T);
  
  int ncorr_2pts=2;
  name corr_name_2pts[2]={"VKVK","P5P5"};
  int npieces_2pts[2]={3,1};

  int ncorr_3pts=20;
  name corr_name_3pts[20]={"S0P5","VKP5","V0P5","TKP5","S0VK","VJVK","VKVJ","P5VK","AKVK","AJVK","AKVJ","A0VK","AKV0","A0V0","TJVK","TKVJ","BKVK","BJVK","BKVJ","BKV0"};

  int npieces_3pts[20]={1,3,1,3,3,3,3,3,3,3,3,3,3,1,3,3,3,3,3,3};
  
  //parse 2pts_30_00
  parse_file("packed_2pts_30_00","2pts_30_00",ncorr_2pts,corr_name_2pts,npieces_2pts);
  
  //parse 2pts_30_30
  parse_file("packed_2pts_30_30","2pts_30_30",ncorr_2pts,corr_name_2pts,npieces_2pts);
  
  //parse 3pts_sp0_30_30
  parse_file("packed_3pts_sp0_30_30","3pts_sp0_30_30",ncorr_3pts,corr_name_3pts,npieces_3pts);
  
  //parse 3pts_sp1_30_30
  parse_file("packed_3pts_sp1_30_30","3pts_sp1_30_30",ncorr_3pts,corr_name_3pts,npieces_3pts);
  
  return 0;
}
