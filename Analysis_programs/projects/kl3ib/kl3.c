#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>
#include "kl3ib_common.cpp"

//////////////////////////////////////////////////////////////////////

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

int main()
{
  //read the input file
  read_input();

  //allocate vector where to load data
  jack_vec *P5_V0_P5[nmass][nmass][nmoms][nmoms];
  jack_vec *P5_P5[nmass][nmoms];
  jack_vec *double_ratio=jack_vec_malloc(T);
  jack_vec *double_ratio_simm=jack_vec_malloc(TH);
  
  //allocate and load two and three point for Zv
  int r=0,im_spec=0;
  for(int im1=0;im1<nmass;im1++)
    for(int ik1=0;ik1<nmoms;ik1++)
      {
	P5_P5[im1][ik1]=jack_vec_malloc(T);
	load_improved_charged_meson_two_points_P5_P5(P5_P5[im1][ik1],base_path,nmoms,nmass,r,im_spec,im1,ik1,theta,L);
	
	for(int im2=0;im2<nmass;im2++)
	  for(int ik2=0;ik2<nmoms;ik2++)
	    {
	      P5_V0_P5[im1][im2][ik1][ik2]=jack_vec_malloc(T);
	      load_improved_charged_meson_three_points_P5_V0_P5(P5_V0_P5[im1][im2][ik1][ik2],base_path,nmoms,nmass,r,im_spec,im1,im2,ik1,ik2,theta,L);
	    }
	}

  //calculate Zv
  jack_vec *two_Zv[nmass];
  for(int imass=0;imass<nmass;imass++)
    {
      two_Zv[imass]=jack_vec_malloc(T);
      jack_frac_jack_vec(two_Zv[imass],P5_P5[imass][0]->data[TH],P5_V0_P5[imass][imass][0][0]);
    }
  
  //fit mass
  jack M[nmass][nmoms],Z;
  for(int imass=0;imass<nmass;imass++) for(int ik=0;ik<nmoms;ik++) jack_fit_mass(Z,M[imass][ik],P5_P5[imass][ik],12,23);
  
  int im_K=0;
  int im_Pi=0;
  
  for(int ik_K=0;ik_K<nmoms;ik_K++)
    for(int ik_Pi=0;ik_Pi<nmoms;ik_Pi++)
      if(theta[ik_K]>=0)
	{
	  //calculate the form factor
	  jack den;
	  jack_prod_jack(den,P5_P5[im_K][ik_K]->data[TH],P5_P5[im_Pi][ik_Pi]->data[TH]);
	  
	  //calculate the num
	  jack_vec_prod_jack_vec(double_ratio,P5_V0_P5[im_K][im_Pi][ik_K][ik_Pi],P5_V0_P5[im_Pi][im_K][ik_Pi][ik_K]);

	  //divide by den
	  jack_vec_fracassign_jack(double_ratio,den);
	  
	  //multiply by 2Zv
	  jack_vec_prodassign_jack_vec(double_ratio,two_Zv[im_K]);
	  jack_vec_prodassign_jack_vec(double_ratio,two_Zv[im_Pi]);
	  
	  //calculate sqrt
	  jack_vec_sqrt_jack_vec(double_ratio,double_ratio);
	  
	  //simmetrize
	  jack_vec_simmetrize(double_ratio_simm,double_ratio,1);
	  
	  //calculate the form factor
	  jack C;
	  constant_fit(C,double_ratio_simm,10,14);
	  //calculate the transferred impulse
	  jack Q2;
	  calculate_Q2(Q2,M[im_K][ik_K],theta[ik_K],M[im_Pi][ik_Pi],theta[ik_Pi],L);
	  
	  printf("%g %g %g %g\n",Q2[njack],C[njack],jack_error(Q2),jack_error(C));
	}
  
  return 0;
}
