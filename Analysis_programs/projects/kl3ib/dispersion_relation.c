#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>
#include "kl3ib_common.cpp"

//////////////////////////////////////////////////////////////////////

int T;
int nmoms,nmass;

void load_averaged_two_points(jack_vec *c,char *base_path,int im1,int im2,int ik1,int ik2,int r1,int r2)
{
  jack_vec *a=jack_vec_malloc(T);
  jack_vec *b=jack_vec_malloc(T);
  
  char path[1024];
  sprintf(path,"%s/oPPo-ss_conf.1.dat",base_path);
  
  int ri=0;
  
  //read the 2 two points with opposite momentum
  read_two_points(a,path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
  read_two_points(b,path,nmoms,nmass,im1,im2,ik2,ik1,r1,r2,ri);
  
  //averages, putting the "-" and the reciprocal spatial volume factor
  jack_vec_summ_jack_vec(c,b,a);
  jack_vec_prodassign_double(c,-0.5/(24*24*24));

  //free vectors
  jack_vec_free(a);
  jack_vec_free(b);
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
  int r1=0,r2=r1;

  //load the pion
  int im1=0,im2=0;
  
  //open all the out file
  FILE *Z_out=open_file("Z_out","w");
  FILE *m_out=open_file("m_out","w");
  FILE *m2_out=open_file("m2_out","w");
  FILE *delta_m2_out=open_file("delta_m2_out","w");
  
  //loop over triangular impulses
  for(int ik1=0;ik1<nmoms;ik1++)
    for(int ik2=ik1;ik2<nmoms;ik2++)
      {
	//calculate q2
	double q2=3*pow((theta[ik1]-theta[ik2])*2*M_PI/24,2);

	//load the correlation functions, averaged	
	load_averaged_two_points(c,base_path,im1,im2,ik1,ik2,r1,r2);
	
	//fit Z and mass
	jack Z,m;
	jack_fit_mass(Z,m,c,12,23);
	
	//calculate squared mass
	jack m2;
	jack_prod_jack(m2,m,m);
	
	//if this is the first mass, copy it to m2_0
	jack m2_0;
	if((ik1==0)&&(ik2==0)) memcpy(m2_0,m2,sizeof(jack));
	
	//subtract m2_0 from m2
	jack delta_m2;
	jack_subt_jack(delta_m2,m2,m2_0);
	
	//print results
	fprintf(Z_out,"%g %g %g\n",q2,Z[njack],jack_error(Z));
	fprintf(m_out,"%g %g %g\n",q2,m[njack],jack_error(m));
	fprintf(m2_out,"%g %g %g\n",q2,m2[njack],jack_error(m2));
	fprintf(delta_m2_out,"%g %g %g\n",q2,delta_m2[njack],jack_error(delta_m2));
      }

  //close output files
  fclose(Z_out);
  fclose(m_out);
  fclose(m2_out);
  fclose(delta_m2_out);
  
  //free the jacknife vector used to load data
  jack_vec_free(c);
  
  return 0;
}
