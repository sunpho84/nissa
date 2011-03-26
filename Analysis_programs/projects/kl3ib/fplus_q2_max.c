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
  jack_vec *P5_Vmu_P5[nmass][nmass][nmoms][nmoms][4];
  jack_vec *P5_P5[nmass][nmoms];
  
  //allocate and load two and three point for Zv
  int r=0,im_spec=0;
  int im_s=1;
  int im_l=0;

  for(int im1=0;im1<nmass;im1++)
    for(int ik1=0;ik1<nmoms;ik1++)
      {
	P5_P5[im1][ik1]=jack_vec_malloc(T);
	load_improved_charged_meson_two_points_P5_P5(P5_P5[im1][ik1],base_path,nmoms,nmass,r,im_spec,im1,ik1,theta,L);
	
	for(int mu=0;mu<4;mu++)
	  for(int im2=0;im2<nmass;im2++)
	    for(int ik2=0;ik2<nmoms;ik2++)
	      {
		P5_Vmu_P5[im1][im2][ik1][ik2][mu]=jack_vec_malloc(T);
		load_improved_charged_meson_three_points_P5_Vmu_P5(P5_Vmu_P5[im1][im2][ik1][ik2][mu],base_path,nmoms,nmass,r,im1,im2,im_spec,ik1,ik2,mu,theta,L);
	      }
	}

  //calculate Zv
  jack_vec *Zv[nmass];
  for(int imass=0;imass<nmass;imass++)
    {
      Zv[imass]=jack_vec_malloc(T);
      jack_frac_jack_vec(Zv[imass],P5_P5[imass][0]->data[TH],P5_Vmu_P5[imass][imass][0][0][0]);
      jack_vec_fracassign_double(Zv[imass],2);
    }
  
  //fit mass
  jack M[nmass][nmoms],Z;
  for(int imass=0;imass<nmass;imass++)
    for(int ik=0;ik<nmoms;ik++)
      jack_mass_fit(Z,M[imass][ik],P5_P5[imass][ik],12,23);
  
  //calculate M2_K-M2_Pi
  jack DM2;
  for(int ij=0;ij<njack+1;ij++)
    DM2[ij]=M[im_s][0][ij]*M[im_s][0][ij]-M[im_l][0][ij]*M[im_l][0][ij];
  
  jack_vec *ratio=jack_vec_malloc(T);
  jack_vec *ratio_simm=jack_vec_malloc(TH);
  
  jack_vec_prod_jack_vec(ratio,P5_Vmu_P5[im_l][im_s][0][0][0],P5_Vmu_P5[im_s][im_l][0][0][0]);
  jack_vec_fracassign_jack_vec(ratio,P5_Vmu_P5[im_s][im_s][0][0][0]);
  jack_vec_fracassign_jack_vec(ratio,P5_Vmu_P5[im_l][im_l][0][0][0]);
  
  //fit f0q2_max
  jack_vec_simmetrize(ratio_simm,ratio,1);
  jack f0q2_max;
  jack_constant_fit(f0q2_max,ratio_simm,10,14);
  printf("%g %g\n",f0q2_max[njack],jack_error(f0q2_max));

  jack_vec_print_to_file("fplus_q2_0",ratio_simm);
  
  return 0;
}
