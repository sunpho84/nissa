#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>
#include "kl3ib_common.cpp"

//////////////////////////////////////////////////////////////////////

int T,TH,L;
int nmoms,nmass;

char base_path[1024];
double *mass,*theta;

void read_input()
{
  FILE *input=open_file("input","r");
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");

  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  L=TH=T/2;

  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list"); 
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");

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

void calculate_d0_A0_P5(jack_vec *d0_A0_P5,jack_vec *A0_P5)
{
  int T=A0_P5->nel;
  jack_vec_check_equal_nel(d0_A0_P5,A0_P5);
  
  for(int ijack=0;ijack<njack+1;ijack++)
    {
      d0_A0_P5->data[0][ijack]=d0_A0_P5->data[T-1][ijack]=0;
      for(int t=1;t<(T-1);t++) d0_A0_P5->data[t][ijack]=(A0_P5->data[t+1][ijack]-A0_P5->data[t-1][ijack])*0.5;
    }
}

int main()
{
  //read the input file
  read_input();

  //allocate vector where to load data
  jack_vec *P5_P5=jack_vec_malloc(T);
  jack_vec *A0_P5=jack_vec_malloc(T);
  jack_vec *d0_A0_P5=jack_vec_malloc(T);
  jack_vec *ratio_averaged=jack_vec_malloc(T);
  jack_vec *ratio_averaged_simm=jack_vec_malloc(TH);
  jack_vec *ratio[2];
  
  //load two points for Zv
  int nr=1;
  for(int r=0;r<nr;r++)
    {
      ratio[r]=jack_vec_malloc(T);

      int im_spec=0,im_valence=0;
      load_standing_charged_meson_two_points_A0_P5(A0_P5,r,im_spec,im_valence);
      load_standing_charged_meson_two_points_P5_P5(P5_P5,r,im_spec,im_valence);
      
      calculate_d0_A0_P5(d0_A0_P5,A0_P5);
      jack_vec_frac_jack_vec(ratio[r],P5_P5,d0_A0_P5);
      jack_vec_prodassign_double(ratio[r],2*mass[0]);
      
      //print results
      if(r==0)
	{
	  print_corr_to_file("P5_P5.out",P5_P5);
	  print_corr_to_file("A0_P5.out",A0_P5);
	  print_corr_to_file("d0_A0_P5.out",d0_A0_P5);
	  print_corr_to_file("ratio.out",ratio[r]);
	}
    }
  
  jack_vec_average(ratio_averaged,ratio,nr);
  jack_vec_simmetrize(ratio_averaged_simm,ratio_averaged,1);
  
  jack Zv;
  constant_fit(Zv,ratio_averaged_simm,6,18);
  printf("Zv: %g %g\n",Zv[njack],jack_error(Zv));
  
  return 0;
}
