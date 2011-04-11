#include <analysis_include.h>
#include "kl3ib_common.cpp"

int T,L,TH;

int ijc;
char base_path[1024];
jack_vec *K_P5P5_simm,*Ks_P5P5_simm,*ratio_P5P5_simm;
jack_vec *K_A0P5_simm,*ratio_A0P5_simm;
double a390=1/2.32;

int nmoms,nmass=3;
double *mass,*theta;
int tmin=12,tmax=23;

double fun_A0P5(double Z_A0P5,double M_A0P5,int t,int TH)
{
  return Z_A0P5*exp(-M_A0P5*TH)*sinh(M_A0P5*(t-TH));
}

double fun_P5P5(double Z_P5P5,double M_P5P5,int t,int TH)
{
  return Z_P5P5*exp(-M_P5P5*TH)*cosh(M_P5P5*(TH-t))/M_P5P5;
}

//return the ratio between the correlation function with and without insertion
double fun_ratio_P5P5(double A,double SL,double M,int t,int TH)
{
  return A+SL*(t-TH)*tanh(M*(t-TH));
}

double fun_ratio_A0P5(double A,double SL,double M,int t,int TH)
{
  return A+SL*(TH-t)/tanh(M*(TH-t));
}

////////////////  fit of slope and ratio of A0P5  /////////////////

//calculate the chi square
double chi2_mass_ratio_A0P5(double A,double SL,double C,double M,int ij)
{
  double ch2=0;
  int TH=K_A0P5_simm->nel;

  for(int t=tmin;t<tmax;t++)
      {
	double ch2_mK=pow((K_A0P5_simm->data[t][ij]-fun_A0P5(C,M,t,TH))/jack_error(K_A0P5_simm->data[t]),2);
	double ch2_ratio=pow((ratio_A0P5_simm->data[t][ij]-fun_ratio_A0P5(A,SL,M,t,TH))/jack_error(ratio_A0P5_simm->data[t]),2);
	
	ch2+=ch2_mK+ch2_ratio;
      }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void ch2_mass_ratio_A0P5_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0];
  double SL=p[1];
  double C=p[2];
  double M=p[3];
  
  ch=chi2_mass_ratio_A0P5(A,SL,C,M,ijc);
}

void jack_fit_mass_and_ratio_A0P5(jack A,jack SL,jack C,jack M)
{
  for(ijc=0;ijc<njack+1;ijc++)
    {
      TMinuit minu;
      
      double meff=-log(K_A0P5_simm->data[tmin+1][ijc]/K_A0P5_simm->data[tmin][ijc]);
      double cestim=fun_A0P5(1,meff,tmin,TH)/K_A0P5_simm->data[tmin][ijc];
      
      minu.SetPrintLevel(-1);
      
      minu.DefineParameter(0,"A",92,0.0001,0,0);
      minu.DefineParameter(1,"SL",1,0.0001,0,0);
      minu.DefineParameter(2,"C",cestim,0.0001,0,0);
      minu.DefineParameter(3,"M",meff,0.0001,0,0);
      
      minu.SetFCN(ch2_mass_ratio_A0P5_wr);
      
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A[ijc],dum);
      minu.GetParameter(1,SL[ijc],dum);
      minu.GetParameter(2,C[ijc],dum);
      minu.GetParameter(3,M[ijc],dum);
    }
}

////////////////  fit of slope and ratio of P5P5  /////////////////

//calculate the chi square
double chi2_mass_ratio_P5P5(double A,double SL,double C,double M,int ij)
{
  double ch2=0;
  int TH=K_P5P5_simm->nel;

  for(int t=tmin;t<tmax;t++)
    {
      double ch2_mK=pow((K_P5P5_simm->data[t][ij]-fun_P5P5(C,M,t,TH))/jack_error(K_P5P5_simm->data[t]),2);
      double ch2_ratio=pow((ratio_P5P5_simm->data[t][ij]-fun_ratio_P5P5(A,SL,M,t,TH))/jack_error(ratio_P5P5_simm->data[t]),2);
      
      ch2+=ch2_mK+ch2_ratio;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void ch2_mass_ratio_P5P5_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0];
  double SL=p[1];
  double C=p[2];
  double M=p[3];
  
  ch=chi2_mass_ratio_P5P5(A,SL,C,M,ijc);
}

void jack_fit_mass_and_ratio_P5P5(jack A,jack SL,jack C,jack M)
{
  for(ijc=0;ijc<njack+1;ijc++)
    {
      TMinuit minu;
      
      double meff=-log(K_P5P5_simm->data[tmin+1][ijc]/K_P5P5_simm->data[tmin][ijc]);
      double cestim=fun_P5P5(1,meff,tmin,TH)/K_P5P5_simm->data[tmin][ijc];
      
      minu.SetPrintLevel(-1);
      
      minu.DefineParameter(0,"A",100,0.0001,0,0);
      minu.DefineParameter(1,"SL",5,0.0001,0,0);
      minu.DefineParameter(2,"C",cestim,0.0001,0,0);
      minu.DefineParameter(3,"M",meff,0.0001,0,0);
      
      minu.SetFCN(ch2_mass_ratio_P5P5_wr);
      
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A[ijc],dum);
      minu.GetParameter(1,SL[ijc],dum);
      minu.GetParameter(2,C[ijc],dum);
      minu.GetParameter(3,M[ijc],dum);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void plot(jack A,jack SL,jack C,jack M)
{
  FILE *fout;
  
  //calculate the correlation function subtracted of the fitted
  jack_vec *K_P5P5_simm_sub=jack_vec_malloc(TH);
  for(int t=0;t<TH;t++)
    for(int ij=0;ij<njack+1;ij++)
      K_P5P5_simm_sub->data[t][ij]=K_P5P5_simm->data[t][ij]-fun_cosh(C[ij],M[ij],t,TH);

  //plot the corr fit
  fout=fopen("fit_corr.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_cosh(C[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,K_P5P5_simm);
  fclose(fout);
  
  //plot the corr fit of Ks
  fout=fopen("fit_corr_ins.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_cosh(C[njack],M[njack],t,TH)/fun_ratio_P5P5(A[njack],SL[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,Ks_P5P5_simm);
  fclose(fout);
  
  //plot the corr subtracted of the fit
  fout=fopen("fit_corr_sub.xmg","w");
  fprintf(fout,"@type xy\n");
  
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,0.0);
  fprintf(fout,"&\n@type xydy\n");
  for(int t=0;t<TH;t++) fprintf(fout,"%d %g %g\n",t,K_P5P5_simm_sub->data[t][njack],jack_error(K_P5P5_simm_sub->data[t]));
  fclose(fout);

  //plot the ratio and the fit
  fout=fopen("fit_ratio.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_ratio_P5P5(A[njack],SL[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,ratio_P5P5_simm);
  fclose(fout);
  
  //plot the ratio minus the fit
  fout=fopen("fit_ratio_sub.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,0.0);
  fprintf(fout,"@type xydy\n");
  jack_vec_fprintf(fout,ratio_P5P5_simm);
  fclose(fout);

  free(K_P5P5_simm_sub);
}

void plot_A0P5(jack A,jack SL,jack C,jack M)
{
  //plot the ratio and the fit
  FILE *fout=fopen("fit_ratio_A0P5.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_ratio_A0P5(A[njack],SL[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,ratio_A0P5_simm);
  fclose(fout);
}

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

int main()
{
  jack_vec *temp_K[2],*temp_Ks[2],*temp_ratio[2];
  
  read_input();
  
  for(int r=0;r<2;r++)
    {
      temp_K[r]=jack_vec_malloc(T);
      temp_Ks[r]=jack_vec_malloc(T);
      temp_ratio[r]=jack_vec_malloc(T);
    }
  
  jack_vec *K_P5P5=jack_vec_malloc(T),*Ks_P5P5=jack_vec_malloc(T);
  jack_vec *K_A0P5=jack_vec_malloc(T),*Ks_A0P5=jack_vec_malloc(T);
  jack_vec *ratio_P5P5=jack_vec_malloc(T),*ratio_A0P5=jack_vec_malloc(T);
  ratio_A0P5_simm=jack_vec_malloc(TH);

  K_P5P5_simm=jack_vec_malloc(TH);
  K_A0P5_simm=jack_vec_malloc(TH);
  Ks_P5P5_simm=jack_vec_malloc(TH);
  ratio_P5P5_simm=jack_vec_malloc(TH);
  
  int im1=0,im2=1; //quello raddoppiato e' il primo!
  int ik1=0,ik2=0;
  int ri=0;
  
  //read P5P5 or A0P5
  for(int ioss=0;ioss<2;ioss++)
    {
      
      //read the 00 and 11 combo
      for(int r=0;r<2;r++)
	{
	  int r1=r,r2=r;
	  
	  char sing[1024],doub[1024];
	  
	  if(ioss==0)
	    {
	      sprintf(sing,"%s/oPPo-ss_conf.1.dat",base_path);
	      sprintf(doub,"%s/oPPo-sd_conf.1.dat",base_path);
	    }
	  else
	    {
	      sprintf(sing,"%s/oA0Po-ss_conf.1.dat",base_path);
	      sprintf(doub,"%s/oA0Po-sd_conf.1.dat",base_path);
	    }
	  
	  read_two_points(temp_K[r],sing,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
	  read_two_points(temp_Ks[r],doub,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
	  
	  jack_vec_prodassign_double(temp_K[r],-1.0/(L*L*L));
	  jack_vec_prodassign_double(temp_Ks[r],-1.0/(L*L*L));
	  
	  //ratio between k with insertion and k
	  jack_vec_frac_jack_vec(temp_ratio[r],temp_Ks[r],temp_K[r]);
	}
      
      //average 00 and 11
      if(ioss==0)
	{
	  jack_vec_average(K_P5P5,temp_K,2);
	  jack_vec_average(Ks_P5P5,temp_Ks,2);
	  jack_vec_average(ratio_P5P5,temp_ratio,2);
	}
      else
	{
	  jack_vec_average(K_A0P5,temp_K,2);
	  jack_vec_average(Ks_A0P5,temp_Ks,2);
	  jack_vec_average(ratio_A0P5,temp_ratio,2);
	}
    }
  
  //simmetrize
  jack_vec_simmetrize(K_P5P5_simm,K_P5P5,1);
  jack_vec_simmetrize(K_A0P5_simm,K_A0P5,-1);
  jack_vec_simmetrize(Ks_P5P5_simm,Ks_P5P5,1);
  jack_vec_simmetrize(ratio_P5P5_simm,ratio_P5P5,1);
  jack_vec_simmetrize(ratio_A0P5_simm,ratio_A0P5,1);
  
  //fit mass and ratio of P5P5
  jack A_P5P5,SL_P5P5,Z_P5P5,M_P5P5;
  jack_fit_mass_and_ratio_P5P5(A_P5P5,SL_P5P5,Z_P5P5,M_P5P5);
  printf("Mass P5P5: %g +- %g = ",M_P5P5[njack],jack_error(M_P5P5));
  printf("%g +- %g GeV\n",M_P5P5[njack]/a390,jack_error(M_P5P5)/a390);
  printf("Slope P5P5: %g %g\n",SL_P5P5[njack],jack_error(SL_P5P5));
  printf("A P5P5: %g %g\n",A_P5P5[njack],jack_error(A_P5P5));
  printf("\n");

  //fit mass and ratio of A0P5
  jack A_A0P5,SL_A0P5,Z_A0P5,M_A0P5;
  jack_fit_mass_and_ratio_A0P5(A_A0P5,SL_A0P5,Z_A0P5,M_A0P5);
  printf("Mass A0P5: %g +- %g = ",M_A0P5[njack],jack_error(M_A0P5));
  printf("%g +- %g GeV\n",M_A0P5[njack]/a390,jack_error(M_A0P5)/a390);
  printf("Slope A0P5: %g %g\n",SL_A0P5[njack],jack_error(SL_A0P5));
  printf("A A0P5: %g %g\n",A_A0P5[njack],jack_error(A_A0P5));
  printf("\n");

  //calculate fK from definition
  jack fK_A0P5;
  for(int ijk=0;ijk<njack+1;ijk++)
    fK_A0P5[ijk]=Z_A0P5[ijk]/sqrt(Z_P5P5[ijk])*0.6108;
  printf("fK (def): %g +- %g = ",fK_A0P5[njack],jack_error(fK_A0P5));
  printf("%g +- %g GeV\n",fK_A0P5[njack]/a390,jack_error(fK_A0P5)/a390);
  
  //calculate fK WI
  jack fK;
  jack_sqrt_jack(fK,Z_P5P5);
  jack_prodassign_double(fK,mass[0]+mass[1]);
  jack_fracassign_jack(fK,M_P5P5);
  jack_fracassign_jack(fK,M_P5P5);
  printf("fK (WI): %g +- %g = ",fK[njack],jack_error(fK));
  printf("%g +- %g GeV\n",fK[njack]/a390,jack_error(fK)/a390);
  printf("\n");

  //calculate delta_fK/fK/delta_m from definition
  jack dfK_fr_fK_2dm;
  for(int ij=0;ij<njack+1;ij++)
    dfK_fr_fK_2dm[ij]=A_A0P5[ij]-0.5*A_P5P5[ij]-0.5*(1/M_P5P5[ij]-TH)*SL_A0P5[ij];
  printf("dfK/fK/2dm (def): %g %g\n",dfK_fr_fK_2dm[njack],jack_error(dfK_fr_fK_2dm));
  
  //calculate delta_fK/fK/delta_m from WI
  jack dfK_fr_fK_2dm_WI;
  for(int ij=0;ij<njack+1;ij++)
    dfK_fr_fK_2dm_WI[ij]=-1/(mass[0]+mass[1])+0.5*(A_P5P5[ij]+(TH-3.0/M_P5P5[ij])*SL_P5P5[ij]);
  printf("dfK/fK/2dm (WI): %g %g\n",dfK_fr_fK_2dm_WI[njack],jack_error(dfK_fr_fK_2dm_WI));
  
  plot(A_P5P5,SL_P5P5,Z_P5P5,M_P5P5);
  plot_A0P5(A_A0P5,SL_A0P5,Z_A0P5,M_A0P5);
  
  return 0;
}
