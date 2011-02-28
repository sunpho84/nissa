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
void load_Breit_degenerate_charged_meson_three_points_P5_V0_P5(jack_vec *P5_V0_P5,int r,int im_spec,int im_valence,int imom)
{
  //load the charged "r" flavour
  int r1=r,r2=!r1,r3=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence,im3=im_valence;
  
  //load the standing meson
  int ik1=imom,ik2=0; //search ik2
  while(theta[ik2]!=-theta[ik1]) ik2++;
  
  //load mu=0
  int mu=0;

  //read
  read_P5_Vmu_P5(P5_V0_P5,base_path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu);

  //Put the 1/spat_vol factor
  jack_vec_prodassign_double(P5_V0_P5,1.0/(L*L*L));
}

//load the mesonic two point function
void load_Breit_charged_meson_two_points_P5_P5(jack_vec *P5_P5,int r,int im_spec,int im_valence,int imom)
{
  //load the charged "r" flavour
  int r1=r,r2=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence;
  
  //load the meson
  int ik1=imom,ik2=0;

  //read
  read_P5_P5(P5_P5,base_path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2);
  
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

int main()
{
  //read the input file
  read_input();

  //allocate vector where to load data
  jack_vec *P5_V0_P5[nmoms];
  jack_vec *P5_P5[nmoms];
  jack_vec *double_ratio=jack_vec_malloc(T);
  jack_vec *double_ratio_simm=jack_vec_malloc(TH);
  
  //allocate and load two and three point for Zv
  int r=0,im_spec=0,im_valence=0;
  for(int imom=0;imom<nmoms;imom++)
    {
      P5_P5[imom]=jack_vec_malloc(T);
      P5_V0_P5[imom]=jack_vec_malloc(T);

      load_Breit_degenerate_charged_meson_three_points_P5_V0_P5(P5_V0_P5[imom],r,im_spec,im_valence,imom);
      load_Breit_charged_meson_two_points_P5_P5(P5_P5[imom],r,im_spec,im_valence,imom);
    }
  
  jack_vec *two_Zv=jack_vec_malloc(T);
  jack_frac_jack_vec(two_Zv,P5_P5[0]->data[TH],P5_V0_P5[0]);  
  
  for(int imom=0;imom<nmoms;imom++)
    if(theta[imom]>=0)
      {
	jack_vec_frac_jack(double_ratio,P5_V0_P5[imom],P5_P5[imom]->data[TH]);
	jack_vec_prodassign_jack_vec(double_ratio,two_Zv);
	jack_vec_simmetrize(double_ratio_simm,double_ratio,1);

	//print the double ratio
	char path_out[1024];
	sprintf(path_out,"EM_ff_mom%d.out",imom);
	print_corr_to_file(path_out,double_ratio_simm);
	
	jack C;
	constant_fit(C,double_ratio_simm,10,14);
	printf("%g %g %g\n",-3*pow(2*2*M_PI/L*theta[imom],2),C[njack],jack_error(C));
      }
  
  //print results
  for(int imom=0;imom<nmoms;imom++)
    {
      char path_out[1024];
      sprintf(path_out,"P5_P5_mom%d.out",imom);
      print_corr_to_file(path_out,P5_P5[imom]);
      sprintf(path_out,"P5_V0_P5_mom%d.out",imom);
      print_corr_to_file(path_out,P5_V0_P5[imom]);
    }

  //free the jacknife vector used to load data
  for(int imom=0;imom<nmoms;imom++)
    {
      jack_vec_free(P5_V0_P5[imom]);
      jack_vec_free(P5_P5[imom]);
    }
  jack_vec_free(double_ratio);
  jack_vec_free(double_ratio_simm);
  jack_vec_free(two_Zv);
  
  return 0;
}
