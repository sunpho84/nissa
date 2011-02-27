#pragma once

//return the ratio between the correlation function with and without insertion
double fun_ratio(double A,double SL,double M,int t,int TH)
{
  return A+SL*(t-TH)*tanh(M*(TH-t));
}

void check_interval(int var,int min,int max)
{
  if(var>max || var<min) crash("Asked for correlation of impossible combination\n",1);
}

//read a particular two point passed as argument
void read_two_points(jack_vec *c,const char *in,int nmoms,int nmass,int im1,int im2,int ik1,int ik2,int r1,int r2,int ri)
{
  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(ri,0,2);
  
  int i=2*(ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2))))+ri;
  jack_vec_load_nazario_format(c,in,i);
}

//read a particular three point passed as argument
void read_three_points(jack_vec *c,const char *in,int nmoms,int nmass,int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu,int ri)
{
  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(im3,0,nmass);check_interval(r3,0,2);
  check_interval(ri,0,2);
  check_interval(mu,0,4);
  
  int i=4*2*(ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2+2*nmass*(im3*2+r3)))))+ri;
  printf("%d\n",i);
  FILE *file=open_file(in,"r");
  int nel=c->nel;

  if(fseek(file,i*sizeof(jack)*nel,SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",i);
      exit(1);
    }

  double data_in[nel][4][2][njack+1];
  int stat=fread(data_in,sizeof(jack)*4*2*nel,1,file);
  if(stat!=1)
    {
      if(stat==EOF) crash("Error, reached EOF while reading data!\n",1);
      else
        {
          perror("Error while reading data");
          exit(1);
        }
    }
  
  for(int i=0;i<nel;i++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[i][ijack]=data_in[i][mu][ri][ijack];
}
