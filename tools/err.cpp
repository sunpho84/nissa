#include <stdio.h>
#include <stdlib.h>

#include "nissa.h"

//count lines
int count_nline(char *path)
{
  FILE *fin=open_file(path,"r");
  int nline=0;
  while(!feof(fin))
    {
      unsigned long len;
      fgetln(fin,&len);
      if(!feof(fin)) nline++;
    }
  fclose(fin);
  
  return nline;
}

//load the file
void load_file(double *in,char *path,int nel_exp)
{
  vector_reset(in);

  FILE *fin=open_file(path,"r");
  
  //read
  for(int iel=0;iel<nel_exp;iel++)
    if(fscanf(fin,"%lg",in+iel)!=1) crash("error reading element %d",iel);
  
  //verify to have reached eof
  //double dum;
  //if(fscanf(fin,"%lg",&dum)==1) crash("error, the file contains more data than lines, read: %lg",dum);
  
  fclose(fin);
}

int main(int narg,char **arg)
{
  init_nissa();
  
  if(narg<2) crash("Use %s file",arg[0]);
  char *path=arg[1];
  
  //count the number of lines
  int nline=count_nline(path);
  master_printf("The file contains %d lines\n",nline);

  //round to next power of 2
  int buf_nel=1;
  while(2*buf_nel<nline) buf_nel*=2;
  master_printf("Buf size: %d\n",buf_nel);
  nline=buf_nel;
  
  //init the grid
  nissa_use_eo_geom=false;
  init_grid(buf_nel,1);
  
  //load the file
  double *in=nissa_malloc("in",nline,double);
  load_file(in,path,nline);
  
  //subtract the average
  double ave=0;
  for(int iel=0;iel<nline;iel++) ave+=in[iel];
  ave/=nline;
  for(int iel=0;iel<nline;iel++) in[iel]-=ave;
  
  //bufferize to compute fourier transform
  complex *buf=nissa_malloc("buf",buf_nel,complex);
  vector_reset(buf);
  for(int iel=0;iel<nline;iel++) buf[iel][RE]=in[iel];
  
  //pass to p space
  fft1d(buf,buf,1,0,1,0);
  
  //square
  for(int iel=0;iel<buf_nel;iel++) safe_complex_conj2_prod(buf[iel],buf[iel],buf[iel]);
  
  //repass to x space
  fft1d(buf,buf,1,0,-1,1);
  
  //unbuffer
  double *out=nissa_malloc("out",nline,double);
  for(int iel=0;iel<nline;iel++) out[iel]=buf[iel][RE];

  //normalize
  double invar=1/out[0];
  for(int iel=0;iel<buf_nel;iel++) out[iel]*=invar;

  
  FILE *fout=open_file("temp","w");
  for(int iel=0;iel<nline;iel++) fprintf(fout,"%16.16lg\n",out[iel]);
  fclose(fout);
  
  nissa_free(out);
  nissa_free(buf);
  nissa_free(in);
  
  close_nissa();
  
  return 0;
}
