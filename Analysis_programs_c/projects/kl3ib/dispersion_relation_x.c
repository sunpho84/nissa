//x^2=p^2/E^2=p^2/(p^2+m^2)->
// -> 1/x^2=1+m^2/p^2 -> m^2=p^2*(1/x^2-1)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>
#include "kl3ib_common.cpp"

//////////////////////////////////////////////////////////////////////

int T,TH,L;
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
  TH=L=T/2;
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  read_formatted_from_file_expecting((char*)&nmoms,input,"%d","nmoms");
  expect_string_from_file(input,"theta_list");
  double *theta=(double*)malloc(sizeof(double)*nmoms);
  for(int imom=0;imom<nmoms;imom++) read_formatted_from_file((char*)&(theta[imom]),input,"%lg","theta");
  fclose(input);

  //allocate vector where to load data
  jack_vec *t0=jack_vec_malloc(T);
  jack_vec *ti=jack_vec_malloc(T);

  //allocate and load two and three point for Zv
  int r=0,im_spec=0;
  int im_s=1;
  int im_l=0;
  
  for(int ik=1;ik<5;ik++)
    {
      //load t1,t2,t3 and average
      load_improved_charged_meson_three_points_P5_Vmu_P5(ti,base_path,nmoms,nmass,r,im_s,im_s,im_spec,ik,ik,1,theta,L);
      load_improved_charged_meson_three_points_P5_Vmu_P5(t0,base_path,nmoms,nmass,r,im_s,im_s,im_spec,ik,ik,2,theta,L);
      jack_vec_summassign_jack_vec(ti,t0);
      load_improved_charged_meson_three_points_P5_Vmu_P5(t0,base_path,nmoms,nmass,r,im_s,im_s,im_spec,ik,ik,3,theta,L);
      jack_vec_summassign_jack_vec(ti,t0);
      jack_vec_fracassign_double(ti,sqrt(3));
      
      //load t0
      load_improved_charged_meson_three_points_P5_Vmu_P5(t0,base_path,nmoms,nmass,r,im_s,im_s,im_spec,ik,ik,0,theta,L);
      
      //calculate x
      jack_vec *x_corr=jack_vec_malloc(T);
      jack_vec *x_corr_simm=jack_vec_malloc(TH);
      jack_vec_frac_jack_vec(x_corr,ti,t0);
      jack_vec_simmetrize(x_corr_simm,x_corr,-1);
      char path[1024];
      sprintf(path,"/tmp/x_%d_%d",im_s,ik);
      jack_vec_print_to_file(path,x_corr_simm);
      
      //fit x
      jack x;
      jack_constant_fit(x,x_corr_simm,10,14);
      //printf("%g %g\n",x[njack],jack_error(x));
      
      //calculate m
      double p=sqrt(3)*2*M_PI*theta[ik]/L,p2=p*p;
      jack m2,e2,m;
      double_frac_jack(e2,1.0,x);     // 1/x = e/p
      jack_frac_jack(e2,e2,x);        // 1/x^2 = (e/p)^2
      jack_prodassign_double(e2,p2);  // (p/x)^2 = e^2
      jack_summ_double(m2,e2,-p2);    // p^2*(1/x^2-1) = m^2
      jack_sqrt_jack(m,m2);            // m
      printf("%g, m = %g +- %g, ",p2,m[njack],jack_error(m));
      printf("e2 = %g +- %g\n",e2[njack],jack_error(e2));
    }
  
  return 0;
}
