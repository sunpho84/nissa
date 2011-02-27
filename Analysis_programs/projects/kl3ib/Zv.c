#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>
#include "kl3ib_common.cpp"

//////////////////////////////////////////////////////////////////////

int T;
int nmoms,nmass;

void read_Vmu(jack_vec *c,const char *base_path,int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu,int ri)
{
  char path[1024];
  sprintf(path,"%s/oPVmuPo-sss_conf.1.dat",base_path);
  
  read_three_points(c,path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,ri);
}

int main()
{
  //read the input file
  FILE *input=open_file("input","r");
  char base_path[1024];
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  read_formatted_from_file_expecting((char*)&nmoms,input,"%d","nmoms");
  expect_string_from_file(input,"theta_list");
  double *theta=(double*)malloc(sizeof(double)*nmoms);
  for(int imom=0;imom<nmoms;imom++) read_formatted_from_file((char*)&(theta[imom]),input,"%lg","theta");
  fclose(input);

  //allocate vector where to load data
  jack_vec *c=jack_vec_malloc(T);

  //load the neutral "up" flavour
  int r1=0,r2=r1,r3=r1;

  //load the pion
  int im1=0,im2=im1,im3=im1;
  
  //load the standing meson
  int ik1=0,ik2=0;
 
  for(int mu=0;mu<4;mu++)
    {
      //open all the out file
      char path_out[1024];
      sprintf(path_out,"three_out_%d",mu);
      FILE *three_out=open_file(path_out,"w");
      
      //load the correlation function
      int ri[4]={0,1,1,1};
      read_Vmu(c,base_path,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,ri[mu]);
      
      //print results
      jack_vec_fprintf(three_out,c);
      
      //close output files
      fclose(three_out);
    }
  
  //free the jacknife vector used to load data
  jack_vec_free(c);
  
  return 0;
}
