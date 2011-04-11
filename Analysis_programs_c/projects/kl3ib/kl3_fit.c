#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>
#include "kl3ib_common.cpp"

//////////////////////////////////////////////////////////////////////

int T,TH,L;
int nmoms,nmass,npos_moms;
int *pos_of_mom,*ind_pos_mom;

char base_path[1024];
double *theta;

int im_spec=0;
int im_s=1,im_l=0;

int _ijack_fit_tot;

jack *Zv;
jack_vec ***P5_P5;
double ***dP5_P5;
jack_vec ******P5_Vmu_P5;
double ******dP5_Vmu_P5;

double fun_P5_P5(double G,double E,int t)
{
  return G*G/(2*E)*(exp(-E*t)+exp(-E*(T-t)));
}  

double fun_P5_Vmu_P5(double G_Pi,double E_Pi,double G_K,double E_K,double F,double Zv_Pi,double Zv_K,int t)
{
  return G_Pi*G_K/(2*E_Pi*2*E_K)/sqrt(Zv_Pi*Zv_K)*F*(exp(-E_Pi*t)*exp(-E_K*(TH-t)));
}  

double chi2_fit_tot(double G_Pi,double G_K,double *E_Pi,double *E_K,double *FQ2_0,double *FQ2_i,int ij)
{
  double y_teo,y_spe,dy,contr_ch2;
  
  double ch2=0;

  //fit P5P5
  for(int imom=0;imom<nmoms;imom++)
    {      
      int ipos=ind_pos_mom[imom];
      
      for(int t=12;t<T-12;t++)
	{
	  //Pi
	  y_teo=fun_P5_P5(G_Pi,E_Pi[ipos],t);
	  y_spe=P5_P5[im_l][imom]->data[t][ij];
	  dy=dP5_P5[im_l][imom][t];
	  
	  contr_ch2=pow((y_teo-y_spe)/dy,2);
	  
	  ch2+=contr_ch2;

	  //K
	  y_teo=fun_P5_P5(G_K,E_K[ipos],t);
	  y_spe=P5_P5[im_s][imom]->data[t][ij];
	  dy=jack_error(P5_P5[im_s][imom]->data[t]);
	  
	  contr_ch2=pow((y_teo-y_spe)/dy,2);
	  
	  ch2+=contr_ch2;
	}
    }

  //Vmu
  for(int mu=0;mu<4;mu++)
    for(int ik_l=0;ik_l<nmoms;ik_l++)
      for(int ik_s=0;ik_s<nmoms;ik_s++)
	if(ik_l!=0 || ik_s!=0)
	{
	  int ip_l=ind_pos_mom[ik_l];
	  int ip_s=ind_pos_mom[ik_s];
	  
	  for(int t=10;t<T;t++)
	    if(t<=14||t>T-14)
	      {
		double F;
		if(mu==0) F=FQ2_0[ik_l*nmoms+ik_s];
		else      F=FQ2_i[ik_l*nmoms+ik_s];
		y_teo=fun_P5_Vmu_P5(G_Pi,G_K,E_Pi[ip_l],E_K[ip_s],F,Zv[im_l][ij],Zv[im_s][ij],t);
		y_spe=P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]->data[t][ij];
		dy=jack_error(P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]->data[t]);
		
		contr_ch2=pow((y_teo-y_spe)/dy,2);
		
		//printf("dy %d %d %d %d %d %d %g\n",im_l,im_s,ik_l,ik_s,mu,t,dy);
		ch2+=contr_ch2;
	      }
	}
  
  //printf("%g\n",ch2);
  
  return ch2;///(nmoms*(T-24)+4*nmoms*nmoms*8);
}

void ch2_fit_tot_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double G_Pi=p[0];
  double G_K=p[1];
  
  ch=chi2_fit_tot(
		  G_Pi,
		  G_K,
		  p+2,
		  p+2+npos_moms,
		  p+2+2*npos_moms,
		  p+2+2*npos_moms+nmoms*nmoms,
		  _ijack_fit_tot);
}

void total_fit(jack *E_Pi,jack *E_K,jack G_Pi,jack G_K,jack **FQ2_0,jack **FQ2_i)
{
  TMinuit minu(100);
  
  minu.DefineParameter(0,"G_Pi",0.1,0.0001,0,0);
  minu.DefineParameter(1,"G_K",0.1,0.0001,0,0);

  char name[1024];
  
  for(int imom=0;imom<npos_moms;imom++)
    {
      sprintf(name,"E_Pi_%d",imom);
      minu.DefineParameter(2+imom,name,0.2,0.0001,0,0);
    }
  for(int imom=0;imom<npos_moms;imom++)
    {
      sprintf(name,"E_K_%d",imom);
      minu.DefineParameter(2+imom+npos_moms,name,0.3,0.0001,0,0);
    }
  
  for(int ik_l=0;ik_l<nmoms;ik_l++)
    for(int ik_s=0;ik_s<nmoms;ik_s++)
      {
	sprintf(name,"FQ2_0_%d_%d",ik_l,ik_s);
	minu.DefineParameter(2+ik_l*nmoms+ik_s+2*npos_moms,name,0.3,0.0001,0,0);
      }
  
  for(int ik_l=0;ik_l<nmoms;ik_l++)
    for(int ik_s=0;ik_s<nmoms;ik_s++)
      {
	sprintf(name,"FQ2_i_%d_%d",ik_l,ik_s);
	minu.DefineParameter(2+nmoms*nmoms+ik_l*nmoms+ik_s+2*npos_moms,name,0.3,0.0001,0,0);
      }
  
  minu.SetFCN(ch2_fit_tot_wr);

  for(_ijack_fit_tot=0;_ijack_fit_tot<njack+1;_ijack_fit_tot++)
    {  
      minu.Migrad();
      double dum;
      minu.GetParameter(0,G_Pi[_ijack_fit_tot],dum);
      minu.GetParameter(1,G_K[_ijack_fit_tot],dum);
      
      for(int ipos=0;ipos<npos_moms;ipos++)
	{
	  minu.GetParameter(ipos+2,E_Pi[ipos][_ijack_fit_tot],dum);
	  minu.GetParameter(ipos+2+npos_moms,E_K[ipos][_ijack_fit_tot],dum);
	}
    }

  FILE *fout;
  for(int imom=0;imom<nmoms;imom++)
    {
      int ipos=ind_pos_mom[imom];
      char path[1024];
      sprintf(path,"/tmp/outPi_%d",imom);
      fout=open_file(path,"w");

      fprintf(fout,"@type xydy\n");
      for(int t=0;t<T;t++)
	fprintf(fout,"%d %g %g\n",t,P5_P5[im_l][imom]->data[t][njack],jack_error(P5_P5[im_l][imom]->data[t]));
      
      fprintf(fout,"@type xy\n");
      for(int t=12;t<T-12;t++)
	fprintf(fout,"%d %g\n",t,fun_P5_P5(G_Pi[njack],E_Pi[ipos][njack],t));
      
      fclose(fout);
    }
}

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

  //for each momentum find that with same module and positive
  pos_of_mom=(int*)malloc(sizeof(int)*nmoms);
  ind_pos_mom=(int*)malloc(sizeof(int)*nmoms);
  npos_moms=0;
  for(int imom=0;imom<nmoms;imom++)
    {
      pos_of_mom[imom]=0;
      while(fabs(theta[pos_of_mom[imom]])!=fabs(theta[imom]) || theta[pos_of_mom[imom]]<0) pos_of_mom[imom]++;

      if(theta[imom]>=0) ind_pos_mom[imom]=npos_moms++;
    }      
  
  //associate negative momenta to positive ones
  for(int imom=0;imom<nmoms;imom++)
    if(theta[imom]<0) 
      ind_pos_mom[imom]=ind_pos_mom[pos_of_mom[imom]];

  //allocate and load two and three point for Zv
  int r=0;

  P5_Vmu_P5=(jack_vec******)malloc(sizeof(jack*****)*nmass);
  dP5_Vmu_P5=(double******)malloc(sizeof(double*****)*nmass);
  for(int im_l=0;im_l<nmass;im_l++)
    {
      P5_Vmu_P5[im_l]=(jack_vec*****)malloc(sizeof(jack_vec****)*nmass);
      dP5_Vmu_P5[im_l]=(double*****)malloc(sizeof(double****)*nmass);
      for(int im_s=0;im_s<nmass;im_s++)
	{
	  P5_Vmu_P5[im_l][im_s]=(jack_vec****)malloc(sizeof(jack_vec***)*nmoms);
	  dP5_Vmu_P5[im_l][im_s]=(double****)malloc(sizeof(double***)*nmoms);
	  for(int ik_l=0;ik_l<nmoms;ik_l++)
	    {
	      P5_Vmu_P5[im_l][im_s][ik_l]=(jack_vec***)malloc(sizeof(jack_vec**)*nmoms);
	      dP5_Vmu_P5[im_l][im_s][ik_l]=(double***)malloc(sizeof(double**)*nmoms);
	      for(int ik_s=0;ik_s<nmoms;ik_s++)
		{
		  P5_Vmu_P5[im_l][im_s][ik_l][ik_s]=(jack_vec**)malloc(sizeof(jack_vec*)*4);
		  dP5_Vmu_P5[im_l][im_s][ik_l][ik_s]=(double**)malloc(sizeof(double*)*4);
		  for(int mu=0;mu<4;mu++)
		    {
		      P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]=jack_vec_malloc(T);
		      dP5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]=(double*)malloc(sizeof(double)*T);
		      for(int t=0;t<T;t++)
			dP5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu][t]=jack_error(P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu]->data[t]);
		    }
		}
	    }
	}
    }
  
  P5_P5=(jack_vec***)malloc(sizeof(jack_vec**)*nmass);
  dP5_P5=(double***)malloc(sizeof(double**)*nmass);
  for(int im1=0;im1<nmass;im1++)
    {
      P5_P5[im1]=(jack_vec**)malloc(sizeof(jack_vec*)*nmoms);
      dP5_P5[im1]=(double**)malloc(sizeof(double*)*nmoms);
      for(int ik1=0;ik1<nmoms;ik1++)
	{
	  P5_P5[im1][ik1]=jack_vec_malloc(T);
	  dP5_P5[im1][ik1]=(double*)malloc(sizeof(double)*T);
	  load_improved_charged_meson_two_points_P5_P5(P5_P5[im1][ik1],base_path,nmoms,nmass,r,im_spec,im1,ik1,theta,L);
	  for(int t=0;t<T;t++) dP5_P5[im1][ik1][t]=jack_error(P5_P5[im1][ik1]->data[t]);
	  for(int mu=0;mu<4;mu++)
	    for(int im2=0;im2<nmass;im2++)
	      for(int ik2=0;ik2<nmoms;ik2++)
		load_improved_charged_meson_three_points_P5_Vmu_P5(P5_Vmu_P5[im1][im2][ik1][ik2][mu],base_path,nmoms,nmass,r,im1,im2,im_spec,ik1,ik2,mu,theta,L);
	}
    }

  //calculate Zv
  jack_vec *Zv_corr[nmass];
  Zv=(jack*)malloc(sizeof(jack)*nmass);
  printf("Zv\n");
  for(int imass=0;imass<nmass;imass++)
    {
      Zv_corr[imass]=jack_vec_malloc(T);
      jack_frac_jack_vec(Zv_corr[imass],P5_P5[imass][0]->data[TH],P5_Vmu_P5[imass][imass][0][0][0]);
      jack_vec_fracassign_double(Zv_corr[imass],2);
      constant_fit(Zv[imass],Zv_corr[imass],6,18);
      
      printf("Zv %d = %g +- %g\n",imass,Zv[imass][njack],jack_error(Zv[imass]));
    }
  
  //Fitted parameters
  jack E_Pi[npos_moms],E_K[npos_moms];
  jack G_Pi,G_K;
  jack *FQ2_0[nmoms],*FQ2_i[nmoms];
  for(int imom=0;imom<nmoms;imom++)
    {
      FQ2_0[imom]=(jack*)malloc(sizeof(jack)*nmoms);
      FQ2_i[imom]=(jack*)malloc(sizeof(jack)*nmoms);
    }

  total_fit(E_Pi,E_K,G_Pi,G_K,FQ2_0,FQ2_i);
  
  return 0;
}
