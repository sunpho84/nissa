#pragma once

#include "include.h"

int T,TH,L;
int nmoms,nmass;

char base_path[1024];
double *theta;

void read_input()
{
  FILE *input=open_file("input","r");
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  L=TH=T/2;
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  read_formatted_from_file_expecting((char*)&nmoms,input,"%d","nmoms");
  expect_string_from_file(input,"theta_list");
  theta=(double*)malloc(sizeof(double)*nmoms);
  for(int imom=0;imom<nmoms;imom++) read_formatted_from_file((char*)&(theta[imom]),input,"%lg","theta");
  fclose(input);
}

void check_interval(int var,int min,int max)
{
  if(var>=max || var<min)
    {
      fprintf(stderr,"Asked for correlation of impossible combination, %d not in the interval: [%d,%d)\n",var,min,max);
      exit(1);
    }
}

//read a particular two point passed as argument
jtvec read_two_points(const char *in,int T,int njack,int nmoms,int nmass,int im1,int im2,int ik1,int ik2,int r1,int r2,int ri)
{
  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(ri,0,2);
  
  int icorr=2*(ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2))))+ri;
  return jvec_load_naz(in,T,njack,icorr);
}
