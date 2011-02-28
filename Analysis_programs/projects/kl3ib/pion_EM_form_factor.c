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

//load the mesonic three point function
void load_improved_degenerate_charged_meson_three_points_P5_V0_P5(jack_vec *P5_V0_P5,int r,int im_spec,int im_valence,int ik1,int ik2)
{
  //load the charged "r" flavour
  int r1=r,r2=!r1,r3=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence,im3=im_valence;
  
  //load mu=0
  int mu=0;

  //read
  read_improved_P5_Vmu_P5(P5_V0_P5,base_path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,theta);

  //Put the 1/spat_vol factor
  jack_vec_prodassign_double(P5_V0_P5,1.0/(L*L*L));
}

//load the mesonic two point function
void load_improved_charged_meson_two_points_P5_P5(jack_vec *P5_P5,int r,int im_spec,int im_valence,int ik1)
{
  //load the charged "r" flavour
  int r1=r,r2=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence;
  
  //load the meson
  int ik2=0;

  //read
  read_improved_P5_P5(P5_P5,base_path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,theta);
  
  //Put the -1/spat_vol factor
  jack_vec_prodassign_double(P5_P5,-1.0/(L*L*L));
}

void print_corr_to_file(const char *path,jack_vec *corr)
{
  //open the out file
  FILE *fout=open_file(path,"w");
  jack_vec_fprintf(fout,corr);  
  fclose(fout);
}

void calculate_Q2(jack Q2,jack E1,double theta1,jack E2,double theta2)
{
  jack DE;
  jack_subt_jack(DE,E1,E2);
  jack_prod_jack(Q2,DE,DE);
  
  double dp=2*M_PI/L*(theta1-theta2);
  jack_subt_double(Q2,Q2,3*dp*dp);
}

int main()
{
  //read the input file
  read_input();

  //allocate vector where to load data
  jack_vec *P5_V0_P5[nmoms][nmoms];
  jack_vec *P5_P5[nmoms];
  jack_vec *double_ratio=jack_vec_malloc(T);
  jack_vec *double_ratio_simm=jack_vec_malloc(TH);
  
  //allocate and load two and three point for Zv
  int r=0,im_spec=0,im_valence=0;
  for(int ik1=0;ik1<nmoms;ik1++)
    {
      P5_P5[ik1]=jack_vec_malloc(T);
      load_improved_charged_meson_two_points_P5_P5(P5_P5[ik1],r,im_spec,im_valence,ik1);

      for(int ik2=0;ik2<nmoms;ik2++)
	{
	  P5_V0_P5[ik1][ik2]=jack_vec_malloc(T);
	  load_improved_degenerate_charged_meson_three_points_P5_V0_P5(P5_V0_P5[ik1][ik2],r,im_spec,im_valence,ik1,ik2);
	}
    }

  //calculate Zv
  jack_vec *two_Zv=jack_vec_malloc(T);
  jack_frac_jack_vec(two_Zv,P5_P5[0]->data[TH],P5_V0_P5[0][0]);
  
  //fit mass
  jack M[nmoms],Z;
  for(int ik=0;ik<nmoms;ik++) jack_fit_mass(Z,M[ik],P5_P5[ik],12,23);
  
  for(int ik1=0;ik1<nmoms;ik1++)
    for(int ik2=0;ik2<nmoms;ik2++)
      if(theta[ik2]>=0)
	{
	  //calculate the form factor
	  jack den;
	  jack_prod_jack(den,P5_P5[ik1]->data[TH],P5_P5[ik2]->data[TH]);
	  jack_sqrt_jack(den,den);
	  jack_vec_frac_jack(double_ratio,P5_V0_P5[ik1][ik2],den);
	  jack_vec_prodassign_jack_vec(double_ratio,two_Zv);
	  jack_vec_simmetrize(double_ratio_simm,double_ratio,1);
	  
	  //print the double ratio
	  char path_out[1024];
	  sprintf(path_out,"EM_ff_mom%d.out",ik1);
	  print_corr_to_file(path_out,double_ratio_simm);
	  
	  //calculate the form factor
	  jack C;
	  constant_fit(C,double_ratio_simm,10,14);
	  //calculate the transferred impulse
	  jack Q2;
	  calculate_Q2(Q2,M[ik1],theta[ik1],M[ik2],theta[ik2]);
	  
	  printf("%g %g %g %g\n",Q2[njack],C[njack],jack_error(Q2),jack_error(C));
	}
  
  //print results
  for(int ik1=0;ik1<nmoms;ik1++)
    {
      int ik2=0;
      while(theta[ik2]!=-theta[ik1]) ik2++;

      char path_out[1024];
      sprintf(path_out,"P5_P5_mom%d.out",ik1);
      print_corr_to_file(path_out,P5_P5[ik1]);
      sprintf(path_out,"P5_V0_P5_mom%d.out",ik1);
      print_corr_to_file(path_out,P5_V0_P5[ik1][ik2]);
    }

  //free the jacknife vector used to load data
  for(int ik1=0;ik1<nmoms;ik1++)
    {
      for(int ik2=0;ik2<nmoms;ik2++) jack_vec_free(P5_V0_P5[ik1][ik2]);
      jack_vec_free(P5_P5[ik1]);
    }
  jack_vec_free(double_ratio);
  jack_vec_free(double_ratio_simm);
  jack_vec_free(two_Zv);
  
  return 0;
}
