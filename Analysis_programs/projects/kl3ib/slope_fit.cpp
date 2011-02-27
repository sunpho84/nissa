#include <analysis_include.h>
#include "kl3ib_common.cpp"

#define T 48

int ijc;
jack_vec *K_simm,*Ks_simm,*ratio_simm;
double a390=1;

int tmin=12,tmax=23;

//calculate the chi square
double chi2(double A,double SL,double C,double M,int ij)
{
  double ch2=0;
  int TH=K_simm->nel;

  for(int t=tmin;t<tmax;t++)
    {
      double ch2_mK=pow((K_simm->data[t][ij]-fun_cosh(C,M,t,TH))/jack_error(K_simm->data[t]),2);
      double ch2_ratio=pow((ratio_simm->data[t][ij]-fun_ratio(A,SL,M,t,TH))/jack_error(ratio_simm->data[t]),2);
      
      ch2+=ch2_mK+ch2_ratio;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void ch2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0];
  double SL=p[1];
  double C=p[2];
  double M=p[3];
  
  ch=chi2(A,SL,C,M,ijc);
}

void jack_fit_mass_and_slope(jack A,jack SL,jack C,jack M)
{
  int TH=K_simm->nel;
 
  for(ijc=0;ijc<njack+1;ijc++)
    {
      TMinuit minu;
      
      double meff=-log(K_simm->data[tmin+1][ijc]/K_simm->data[tmin][ijc]);
      double cestim=fun_cosh(1,meff,tmin,TH)/K_simm->data[tmin][ijc];
      
      minu.DefineParameter(0,"A",100,0.0001,0,0);
      minu.DefineParameter(1,"SL",5,0.0001,0,0);
      minu.DefineParameter(2,"C",cestim,0.0001,0,0);
      minu.DefineParameter(3,"M",meff,0.0001,0,0);
      
      minu.SetFCN(ch2_wr);
      
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A[ijc],dum);
      minu.GetParameter(1,SL[ijc],dum);
      minu.GetParameter(2,C[ijc],dum);
      minu.GetParameter(3,M[ijc],dum);
    }
}

void plot(jack A,jack SL,jack C,jack M)
{
  FILE *fout;
  
  int TH=T/2;
  
  //calculate the correlation function subtracted of the fitted
  jack_vec *K_simm_sub=jack_vec_malloc(TH);
  for(int t=0;t<TH;t++)
    for(int ij=0;ij<njack+1;ij++)
      K_simm_sub->data[t][ij]=K_simm->data[t][ij]-fun_cosh(C[ij],M[ij],t,TH);

  //plot the corr fit
  fout=fopen("fit_corr.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_cosh(C[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,K_simm);
  fclose(fout);
  
  //plot the corr fit of Ks
  fout=fopen("fit_corr_ins.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_cosh(C[njack],M[njack],t,TH)/fun_ratio(A[njack],SL[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,Ks_simm);
  fclose(fout);
  
  //plot the corr subtracted of the fit
  fout=fopen("fit_corr_sub.xmg","w");
  fprintf(fout,"@type xy\n");
  
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,0.0);
  fprintf(fout,"&\n@type xydy\n");
  for(int t=0;t<TH;t++) fprintf(fout,"%d %g %g\n",t,K_simm_sub->data[t][njack],jack_error(K_simm_sub->data[t]));
  fclose(fout);

  //plot the ratio and the fit
  fout=fopen("fit_ratio.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,fun_ratio(A[njack],SL[njack],M[njack],t,TH));
  fprintf(fout,"&\n@type xydy\n");
  jack_vec_fprintf(fout,ratio_simm);
  fclose(fout);
  
  //plot the ratio minus the fit
  fout=fopen("fit_ratio_sub.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,0.0);
  fprintf(fout,"@type xydy\n");
  jack_vec_fprintf(fout,ratio_simm);
  fclose(fout);

  /*
  //plot the effective mass
  fout=fopen("effective_mass.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,M[njack]);
  fprintf(fout,"@type xydy\n");
  for(int t=0;t<TH-1;t++) fprintf(fout,"%d %g %g\n",t,eff_mass[t][njack],jack_error(eff_mass[t]));
  fclose(fout);

  //plot the effective mass subtracted
  fout=fopen("effective_mass_sub.xmg","w");
  fprintf(fout,"@type xy\n");
  for(int t=tmin;t<tmax;t++) fprintf(fout,"%d %g\n",t,0.0);
  fprintf(fout,"@type xydy\n");
  for(int t=0;t<TH-1;t++) fprintf(fout,"%d %g %g\n",t,eff_mass_sub[t][njack],jack_error(eff_mass_sub[t]));
  fclose(fout);
  */
  free(K_simm_sub);
}

int main()
{
  int nmoms=5,nmass=3;
  jack_vec *temp_K[2],*temp_Ks[2],*temp_ratio[2];
  
  for(int r=0;r<2;r++)
    {
      temp_K[r]=jack_vec_malloc(T);
      temp_Ks[r]=jack_vec_malloc(T);
      temp_ratio[r]=jack_vec_malloc(T);
    }
  
  jack_vec *K=jack_vec_malloc(T),*Ks=jack_vec_malloc(T),*ratio=jack_vec_malloc(T);
  
  int im1=0,im2=1; //quello raddoppiato e' il primo!
  int ik1=0,ik2=0;
  int ri=0;
  int TH=T/2;

  //read the 00 and 11 combo
  for(int r=0;r<2;r++)
    {
      int r1=r,r2=r;
      
      read_two_points(temp_K[r],"oPPo-ss_conf.1.dat",nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
      read_two_points(temp_Ks[r],"oPPo-sd_conf.1.dat",nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
      
      jack_vec_prodassign_double(temp_K[r],-1.0/(T*TH*TH));
      jack_vec_prodassign_double(temp_Ks[r],-1.0/(TH*TH*TH));
      
      //ratio between k with insertion and k
      jack_vec_frac_jack_vec(temp_ratio[r],temp_Ks[r],temp_K[r]);
    }

  //average 00 and 11
  jack_vec_average(K,temp_K,2);
  jack_vec_average(Ks,temp_Ks,2);
  jack_vec_average(ratio,temp_ratio,2);

  //average left and right
  K_simm=jack_vec_malloc(TH);
  Ks_simm=jack_vec_malloc(TH);
  ratio_simm=jack_vec_malloc(TH);

  jack_vec_simmetrize(K_simm,K,1);
  jack_vec_simmetrize(Ks_simm,Ks,1);
  jack_vec_simmetrize(ratio_simm,ratio,1);
  
  jack A,SL,C,M;
  jack_fit_mass_and_slope(A,SL,C,M);
  jack_vec_fprintf(stdout,K_simm);
  
  printf("Mass: %g %g\n",M[njack]/a390,jack_error(M)/a390);
  printf("Slope: %g %g\n",SL[njack],jack_error(SL));
  
  plot(A,SL,C,M);
  
  return 0;
}
