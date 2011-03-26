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
  jack_vec *double_ratio_simm[2];
  jack_vec *temp=jack_vec_malloc(T);
  for(int mu=0;mu<2;mu++) double_ratio_simm[mu]=jack_vec_malloc(TH);
  
  //allocate and load two and three point for Zv
  int im_spec=0;
  int im_s=1;
  int im_l=0;

  for(int im1=0;im1<nmass;im1++)
    for(int ik1=0;ik1<nmoms;ik1++)
      {
	P5_P5[im1][ik1]=jack_vec_malloc(T);
	
	//average r 0 and 1
	load_improved_charged_meson_two_points_P5_P5(P5_P5[im1][ik1],base_path,nmoms,nmass,0,im_spec,im1,ik1,theta,L);
	load_improved_charged_meson_two_points_P5_P5(temp,base_path,nmoms,nmass,1,im_spec,im1,ik1,theta,L);
	jack_vec_summ_jack_vec(P5_P5[im1][ik1],P5_P5[im1][ik1],temp);
	jack_vec_prod_double(P5_P5[im1][ik1],P5_P5[im1][ik1],0.5);

	for(int mu=0;mu<4;mu++)
	  for(int im2=0;im2<nmass;im2++)
	    for(int ik2=0;ik2<nmoms;ik2++)
	      {
		P5_Vmu_P5[im1][im2][ik1][ik2][mu]=jack_vec_malloc(T);
		//average r = 0 and 1
		load_improved_charged_meson_three_points_P5_Vmu_P5(P5_Vmu_P5[im1][im2][ik1][ik2][mu],base_path,nmoms,nmass,0,im1,im2,im_spec,ik1,ik2,mu,theta,L);
		load_improved_charged_meson_three_points_P5_Vmu_P5(temp,base_path,nmoms,nmass,1,im1,im2,im_spec,ik1,ik2,mu,theta,L);
		jack_vec_summ_jack_vec(P5_Vmu_P5[im1][im2][ik1][ik2][mu],P5_Vmu_P5[im1][im2][ik1][ik2][mu],temp);
		//average k-pi with pi-k
		load_improved_charged_meson_three_points_P5_Vmu_P5(temp,base_path,nmoms,nmass,0,im2,im1,im_spec,ik2,ik1,mu,theta,L);
		jack_vec_summ_jack_vec(P5_Vmu_P5[im1][im2][ik1][ik2][mu],P5_Vmu_P5[im1][im2][ik1][ik2][mu],temp);
		load_improved_charged_meson_three_points_P5_Vmu_P5(temp,base_path,nmoms,nmass,1,im2,im1,im_spec,ik2,ik1,mu,theta,L);
		jack_vec_summ_jack_vec(P5_Vmu_P5[im1][im2][ik1][ik2][mu],P5_Vmu_P5[im1][im2][ik1][ik2][mu],temp);
		jack_vec_prod_double(P5_Vmu_P5[im1][im2][ik1][ik2][mu],P5_Vmu_P5[im1][im2][ik1][ik2][mu],0.25);
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
  jack E[nmass][nmoms],Z;
  for(int imass=0;imass<nmass;imass++)
    {
      int ik=0;
      jack_mass_fit(Z,E[imass][ik],P5_P5[imass][ik],12,23);
      for(int ik=1;ik<nmoms;ik++)
	{
	  double p=2*M_PI*theta[ik]/L;
	  for(int ij=0;ij<njack+1;ij++) E[imass][ik][ij]=sqrt(pow(E[imass][0][ij],2)+p*p*3);
	}
    }

  //calculate M2_K-M2_Pi
  jack DM2;
  for(int ij=0;ij<njack+1;ij++)
    DM2[ij]=E[im_s][0][ij]*E[im_s][0][ij]-E[im_l][0][ij]*E[im_l][0][ij];
  
  //Form factor
  FILE *file_f0=open_file("f0","w");
  FILE *file_fp=open_file("fp","w");
  FILE *file_fm=open_file("fm","w");
  
  jack_vec *temp1=jack_vec_malloc(TH);
  jack_vec *temp2=jack_vec_malloc(TH);

  jack V0,Vi;
  for(int ik_l=0;ik_l<nmoms;ik_l++)
    for(int ik_s=0;ik_s<nmoms;ik_s++)
      {
	printf("-------\n");

	for(int t=0;t<TH;t++)
	  for(int ij=0;ij<=njack;ij++)
	    {
	      double a,b;

	      //temporal direction: simmetrize t with T-t
	      a=(P5_Vmu_P5[im_s][im_l][ik_s][ik_l][0]->data[t][ij]-
		 P5_Vmu_P5[im_s][im_l][ik_s][ik_l][0]->data[T-t][ij])/2;
	      b=(P5_Vmu_P5[im_l][im_s][ik_l][ik_s][0]->data[t][ij]-
		 P5_Vmu_P5[im_l][im_s][ik_l][ik_s][0]->data[T-t][ij])/2;
	      double_ratio_simm[0]->data[t][ij]=a*b;
	
	      //spatial direction of K->Pi
	      a=b=0;
	      for(int mu=1;mu<4;mu++)
		{
		  a+=P5_Vmu_P5[im_s][im_l][ik_s][ik_l][mu]->data[t][ij]+
		    P5_Vmu_P5[im_s][im_l][ik_s][ik_l][mu]->data[T-t][ij];
		  b+=P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]->data[t][ij]+
		    P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]->data[T-t][ij];
		}
	      a/=6;
	      b/=6;
	      
	      //product of K->Pi Pi->K
	      double_ratio_simm[1]->data[t][ij]=a*b;
	    }
	
	//calculate the denominator
	jack_vec_simmetrize(temp1,P5_Vmu_P5[im_l][im_l][ik_l][ik_l][0],-1);
	jack_vec_simmetrize(temp2,P5_Vmu_P5[im_s][im_s][ik_s][ik_s][0],-1);
	jack_vec_prodassign_jack_vec(temp1,temp2);

	//calculate the two ratios
	jack sq;
	jack_prod_jack(sq,E[im_s][ik_s],E[im_l][ik_l]);
	jack_sqrt_jack(sq,sq);
	jack_prodassign_double(sq,2);
	for(int mu=0;mu<2;mu++)
	  {
	    jack_vec_fracassign_jack_vec(double_ratio_simm[mu],temp1);
	    jack_vec_sqrt_jack_vec(double_ratio_simm[mu],double_ratio_simm[mu]);
	    
	    //reassign the correct sign
	    if(P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]->data[12][njack]<0)
	      jack_vec_prodassign_double(double_ratio_simm[mu],-1);

	    //multiply by 2*sqrt(E_K*E_Pi)
	    jack_vec_prodassign_jack(double_ratio_simm[mu],sq);
	  }
	
	//calculate the form factor
	jack_constant_fit(V0,double_ratio_simm[0],10,14);
	jack_constant_fit(Vi,double_ratio_simm[1],10,14);

	printf("V0= %g %g\n",V0[njack]/sq[njack],jack_error(V0)/sq[njack]);
	printf("Vi= %g %g\n",Vi[njack]/sq[njack],jack_error(Vi)/sq[njack]);
	
	char path[1024];
	FILE *fout;
	for(int mu=0;mu<2;mu++)
	  {
	    sprintf(path,"out_double_ratio_simm_P5_V%d_P5_%d_%d",mu,ik_s,ik_l);
	    fout=open_file(path,"w");
	    jack_vec_fprintf(fout,double_ratio_simm[mu]);
	    fclose(fout);
	  }

	//calculate the impulses
	jack Q[4],P[4];
	jack Q2;
	calculate_Q_P(Q,P,E[im_s][ik_s],theta[ik_s],E[im_l][ik_l],theta[ik_l],L);
	quad_jack_prod_quad_jack(Q2,Q,Q);

	if(theta[ik_s]!=-theta[ik_l] && ik_s!=ik_l)
	  {
	    //calculate f+,fm,f0
	    jack fp,fm,f0;
	    for(int ij=0;ij<njack+1;ij++)
	      {
		double P0=P[0][ij],Pi=P[1][ij];
		double Q0=Q[0][ij],Qi=Q[1][ij];
		double DET=P0*Qi-Q0*Pi;
		
		if(ij==njack) printf("%d %d Q0=%g Qi=%g, P0=%g Pi=%g\n",ik_s,ik_l,Q0,Qi,P0,Pi);
		
		//calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
		fp[ij]=(V0[ij]*Qi-Vi[ij]*Q0)/DET;
		fm[ij]=(P0*Vi[ij]-Pi*V0[ij])/DET;
		
		f0[ij]=fp[ij]+Q2[ij]/DM2[ij]*fm[ij];
	      }

	    	    
	    fprintf(file_fp,"%g %g %g %d %d\n",Q2[njack],fp[njack],jack_error(fp),ik_s,ik_l);
	    fprintf(file_fm,"%g %g %g %d %d\n",Q2[njack],fm[njack],jack_error(fm),ik_s,ik_l);
	    fprintf(file_f0,"%g %g %g %d %d\n",Q2[njack],f0[njack],jack_error(f0),ik_s,ik_l);
	  }

	if(theta[ik_s]==-theta[ik_l] && ik_s!=ik_l)
	  {
	    //calculate f+,fm,f0
	    jack fp,fm,f0;
	    for(int ij=0;ij<njack+1;ij++)
	      {
		double P0=P[0][ij],Pi=P[1][ij];
		double Q0=Q[0][ij],Qi=Q[1][ij];
		
		if(ij==njack) printf("%d %d Q0=%g Qi=%g, P0=%g Pi=%g\n",ik_s,ik_l,Q0,Qi,P0,Pi);
		
		//calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
		fm[ij]=Vi[ij]/Qi;
		fp[ij]=(V0[ij]-Q0*fm[ij])/P0;
		
		f0[ij]=fp[ij]+Q2[ij]/DM2[ij]*fm[ij];
	      }
	    	    
	    fprintf(file_fp,"%g %g %g %d %d b\n",Q2[njack],fp[njack],jack_error(fp),ik_s,ik_l);
	    fprintf(file_fm,"%g %g %g %d %d b\n",Q2[njack],fm[njack],jack_error(fm),ik_s,ik_l);
	    fprintf(file_f0,"%g %g %g %d %d\n",Q2[njack],f0[njack],jack_error(f0),ik_s,ik_l);
	  }
	}
  
  return 0;
}
