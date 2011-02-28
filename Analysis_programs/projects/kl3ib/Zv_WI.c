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

//load the mesonic two point function for the Zv
void load_standing_charged_meson_two_points_P5_P5(jack_vec *P5_P5,int r,int im_spec,int im_valence)
{
  //load the charged "r" flavour
  int r1=r,r2=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence;
  
  //load the standing meson
  int ik1=0,ik2=0;
 
  //read
  read_P5_P5(P5_P5,base_path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2);
  
  //Put the -1/spat_vol factor
  jack_vec_prodassign_double(P5_P5,-1.0/(L*L*L));
}

void load_standing_charged_meson_two_points_A0_P5(jack_vec *A0_P5,int r,int im_spec,int im_valence)
{
  //load the charged "r" flavour
  int r1=r,r2=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence;
  
  //load the standing meson
  int ik1=0,ik2=0;
 
  //read
  read_A0_P5(A0_P5,base_path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2);
  
  //Put the -1/spat_vol factor
  jack_vec_prodassign_double(A0_P5,-1.0/(L*L*L));
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
  jack_vec *P5_P5=jack_vec_malloc(T);
  jack_vec *A0_P5=jack_vec_malloc(T);
  jack_vec *ratio=jack_vec_malloc(T);
  jack_vec *ratio_simm=jack_vec_malloc(TH);
  
  //load two points for Zv
  int r=1,im_spec=0,im_valence=0;
  load_standing_charged_meson_two_points_A0_P5(A0_P5,r,im_spec,im_valence);
  load_standing_charged_meson_two_points_P5_P5(P5_P5,r,im_spec,im_valence);
  
  jack_vec_frac_jack_vec(ratio,P5_P5,A0_P5);
  jack_vec_prodassign_double(ratio,-2*0.0064);
  jack_vec_simmetrize(ratio_simm,ratio,-1);
  
  //print results
  print_corr_to_file("P5_P5_out",P5_P5);
  print_corr_to_file("A0_P5_out",A0_P5);
  print_corr_to_file("ratio_out",ratio);
  print_corr_to_file("ratio_simm_out",ratio_simm);
  
  //free the jacknife vector used to load data
  jack_vec_free(P5_P5);
  jack_vec_free(A0_P5);
  jack_vec_free(ratio);
  jack_vec_free(ratio_simm);
  
  return 0;
}
