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

//load the mesonic three point function for the Zv
void load_standing_degenerate_charged_meson_three_points_P5_V0_P5(jack_vec *P5_V0_P5,int r,int im_spec,int im_valence)
{
  //load the charged "r" flavour
  int r1=r,r2=!r1,r3=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence,im3=im_valence;
  
  //load the standing meson
  int ik1=0,ik2=0;
 
  //load mu=0
  int mu=0;

  //read
  read_P5_Vmu_P5(P5_V0_P5,base_path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu);

  //Put the 1/spat_vol factor
  jack_vec_prodassign_double(P5_V0_P5,1.0/(L*L*L));
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
  jack_vec *P5_V0_P5=jack_vec_malloc(T);
  jack_vec *P5_V0_P5_simm=jack_vec_malloc(TH);
  jack_vec *P5_P5=jack_vec_malloc(T);
  jack_vec *Zv_corr=jack_vec_malloc(TH);

  //load two and three point for Zv
  int r=0,im_spec=0,im_valence=0;
  load_standing_degenerate_charged_meson_three_points_P5_V0_P5(P5_V0_P5,r,im_spec,im_valence);
  load_standing_charged_meson_two_points_P5_P5(P5_P5,r,im_spec,im_valence);

  //calculate Zv_corr
  jack_vec_simmetrize(P5_V0_P5_simm,P5_V0_P5,-1);
  jack_frac_jack_vec(Zv_corr,P5_P5->data[TH],P5_V0_P5_simm);
  jack_vec_prodassign_double(Zv_corr,0.5);
  
  //fit constant over Zv
  jack Zv_const;
  constant_fit(Zv_const,Zv_corr,6,18);
  printf("Zv: %g %g\n",Zv_const[njack],jack_error(Zv_const));

  //print results
  print_corr_to_file("P5_P5_out",P5_P5);
  print_corr_to_file("P5_V0_P5_out",P5_V0_P5);
  print_corr_to_file("P5_P5_simm_out",P5_V0_P5_simm);
  print_corr_to_file("Zv_corr_out",Zv_corr);

  //free the jacknife vector used to load data
  jack_vec_free(P5_V0_P5_simm);
  jack_vec_free(P5_V0_P5);
  jack_vec_free(P5_P5);
  jack_vec_free(Zv_corr);
  
  return 0;
}
