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
  jack_vec *double_ratio=jack_vec_malloc(T);
  jack_vec *double_ratio_simm=jack_vec_malloc(TH);
  
  //allocate and load two and three point for Zv
  int r=0,im_spec=0;
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
		load_improved_charged_meson_three_points_P5_Vmu_P5(P5_Vmu_P5[im1][im2][ik1][ik2][mu],base_path,nmoms,nmass,r,im_spec,im1,im2,ik1,ik2,mu,theta,L);
	      }
	}

  //calculate Zv
  jack_vec *two_Zv[nmass];
  for(int imass=0;imass<nmass;imass++)
    {
      two_Zv[imass]=jack_vec_malloc(T);
      jack_frac_jack_vec(two_Zv[imass],P5_P5[imass][0]->data[TH],P5_Vmu_P5[imass][imass][0][0][0]);
    }
  
  //fit mass
  jack M[nmass][nmoms],Z;
  for(int imass=0;imass<nmass;imass++) for(int ik=0;ik<nmoms;ik++) jack_fit_mass(Z,M[imass][ik],P5_P5[imass][ik],12,23);
  
  int im_K=1;
  int im_Pi=0;
  
  //calculate M_K^2-M_Pi^2
  jack squared_mass_diff;
  jack M2_K,M2_Pi;
  jack_prod_jack(M2_K,M[im_K][0],M[im_K][0]);
  jack_prod_jack(M2_Pi,M[im_Pi][0],M[im_Pi][0]);
  jack_subt_jack(squared_mass_diff,M2_K,M2_Pi);
  printf("Squared mass: %g\n",squared_mass_diff[njack]);
  
  //Form factor
  jack FQ2[nmoms][nmoms][4];
  for(int ik_K=0;ik_K<nmoms;ik_K++)
    for(int ik_Pi=0;ik_Pi<nmoms;ik_Pi++)
      if(theta[ik_K]>=0)
	{
	  for(int mu=0;mu<4;mu++)
	    {
	      //calculate the form factor
	      jack den;
	      jack_prod_jack(den,P5_P5[im_K][ik_K]->data[TH],P5_P5[im_Pi][ik_Pi]->data[TH]);
	      
	      //calculate the num
	      jack_vec_prod_jack_vec(double_ratio,P5_Vmu_P5[im_K][im_Pi][ik_K][ik_Pi][mu],P5_Vmu_P5[im_Pi][im_K][ik_Pi][ik_K][mu]);
	      
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
	      constant_fit(FQ2[ik_Pi][ik_K][mu],double_ratio_simm,10,14);
	    }
	  
	  //calculate the impulses
	  jack Q[4],P[4],P2,Q2;
	  calculate_Q_P(Q,P,M[im_K][ik_K],theta[ik_K],M[im_Pi][ik_Pi],theta[ik_Pi],L);
	  quad_jack_prod_quad_jack(P2,P,P);
	  quad_jack_prod_quad_jack(Q2,Q,Q);
	  
	  //calculate (f+)*P2 and (f-)*Q2
	  jack fp_prod_P2,fm_prod_Q2;
	  quad_jack_prod_quad_jack(fp_prod_P2,FQ2[ik_Pi][ik_K],P);
	  quad_jack_prod_quad_jack(fm_prod_Q2,FQ2[ik_Pi][ik_K],Q);
	  //calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
	  jack fp,fm_prime;
	  jack_frac_jack(fp,fp_prod_P2,P2);
	  jack_frac_jack(fm_prime,fm_prod_Q2,squared_mass_diff);
	  
	  //calculate f0
	  jack f0;
	  jack_summ_jack(f0,fp,fm_prime);
	  
	  //calculate fp*P[0]
	  jack fp_prod_P0;
	  jack_prod_jack(fp_prod_P0,fp,P[0]);
	  
	  //if(P[1][njack]==0)
	    {
	      /*printf("P2: %g\n",P2[njack]); 
	      printf("fp_prod_P2: %g\n",fp_prod_P2[njack]);
	      printf("FQ2[0]: %g\n",FQ2[ik_Pi][ik_K][0][njack]);
	      printf("fp: %g\n",fp[njack]); 
	      printf("P[0]: %g\n",P[0][njack]); */
	      printf("%g %g %g %g\n",Q2[njack],f0[njack],jack_error(Q2),jack_error(f0));
	    }
	}
  
  return 0;
}
