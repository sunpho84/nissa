#include "include.h"
#include "kl3_common.cpp"

double fun_A0P5(double Z_A0P5,double M_A0P5,int t){return Z_A0P5*exp(-M_A0P5*TH)*sinh(M_A0P5*(t-TH));}
double fun_P5P5(double Z_P5P5,double M_P5P5,int t){return Z_P5P5*exp(-M_P5P5*TH)*cosh(M_P5P5*(TH-t))/M_P5P5;}
double fun_ratio_P5P5(double A,double SL,double M,int t){return A+SL*(t-TH)*tanh(M*(t-TH));}
double fun_ratio_A0P5(double A,double SL,double M,int t){return A+SL*(TH-t)/tanh(M*(TH-t));}

/////////////////  fit of slope and ratio of A0P5  /////////////////

double (*fit_fun_K)(double,double,int),(*fit_fun_ratio)(double,double,double,int);
int tmin_fit,tmax_fit;
double *K_fit,*ratio_fit;
double *dK_fit,*dratio_fit;
double a390=1/2.32/1000;
double Zv390=0.6108;
double Zp390=0.437;

double dM_K_fis=-6;
double ml=3.6;

//calculate the chi square
double chi2_mass_ratio(double A,double SL,double C,double M)
{
  double ch2=0;

  for(int t=tmin_fit;t<=tmax_fit;t++)
    {
      double ch2_mK=pow((K_fit[t]-fit_fun_K(C,M,t))/dK_fit[t],2);
      double ch2_ratio=pow((ratio_fit[t]-fit_fun_ratio(A,SL,M,t))/dratio_fit[t],2);
      ch2+=ch2_mK+ch2_ratio;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void ch2_mass_ratio_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0];double SL=p[1];double C=p[2];double M=p[3];
  ch=chi2_mass_ratio(A,SL,C,M);
}

void jack_fit_mass_and_ratio(jack &A,jack &SL,jack &C,jack &M,jvec &K,jvec &ratio,int tmin,int tmax,double (*fun_K)(double,double,int),double (*fun_ratio)(double,double,double,int))
{
  if(tmin<0||tmax>TH)
    {
      cerr<<"Error, requiring interval "<<tmin<<" "<<tmax<<endl;
      exit(1);
    }
  tmin_fit=tmin;
  tmax_fit=tmax;

  fit_fun_K=fun_K;
  fit_fun_ratio=fun_ratio;
  
  K_fit=new double[TH+1];
  dK_fit=new double[TH+1];
  ratio_fit=new double[TH+1];
  dratio_fit=new double[TH+1];
  
  //calculate errors
  for(int t=tmin;t<=tmax;t++)
    {
      dK_fit[t]=K.data[t].err();
      dratio_fit[t]=ratio.data[t].err(); 
    }
  
  for(int ijack=0;ijack<njack+1;ijack++)
    {
      //copy data so that glob function may access it
      for(int t=tmin;t<=tmax;t++)
	{
	  K_fit[t]=K.data[t].data[ijack];
	  ratio_fit[t]=ratio.data[t].data[ijack];
	}

      TMinuit minu;
      
      double meff=-log(K_fit[tmin+1]/K_fit[tmin]);
      double cestim=fun_K(1,meff,tmin)/K_fit[tmin];
      
      minu.SetPrintLevel(-1);
      
      minu.DefineParameter(0,"A",100,0.0001,0,0);
      minu.DefineParameter(1,"SL",5,0.0001,0,0);
      minu.DefineParameter(2,"C",cestim,0.0001,0,0);
      minu.DefineParameter(3,"M",meff,0.0001,0,0);
      
      minu.SetFCN(ch2_mass_ratio_wr);
      
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A.data[ijack],dum);
      minu.GetParameter(1,SL.data[ijack],dum);
      minu.GetParameter(2,C.data[ijack],dum);
      minu.GetParameter(3,M.data[ijack],dum);
    }
}

void jack_fit_mass_and_ratio_P5P5(jack &A,jack &SL,jack &C,jack &M,jvec &K,jvec &ratio,int tmin,int tmax)
{jack_fit_mass_and_ratio(A,SL,C,M,K,ratio,tmin,tmax,fun_P5P5,fun_ratio_P5P5);}

void jack_fit_mass_and_ratio_A0P5(jack &A,jack &SL,jack &C,jack &M,jvec &K,jvec &ratio,int tmin,int tmax)
{jack_fit_mass_and_ratio(A,SL,C,M,K,ratio,tmin,tmax,fun_A0P5,fun_ratio_A0P5);}


jvec load_chaveraged_2pts(const char *nome_file,int im1,int im2,int ik1,int ik2,int ri)
{
  char path[1024];
  sprintf(path,"%s/%s",base_path,nome_file);  
  
  //read the 2 two points with opposite r
  return (read_two_points(path,im1,im2,ik1,ik2,0,0,ri)+
	  read_two_points(path,im1,im2,ik1,ik2,1,1,ri))/2;
}

int main()
{
  njack=10;
  read_input();
  
  int im1=0,im2=1; //quello raddoppiato e' il primo!
  int ik1=0,ik2=0;
  int ri=0;
  
  //load the data
  jvec Pi_P5P5 =load_chaveraged_2pts("oPPo-ss_conf.1.dat",im1,im1,ik1,ik2,ri).simmetrized(1);
  jvec K_A0P5 =load_chaveraged_2pts("oA0Po-ss_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(-1);
  jvec Ks_A0P5=load_chaveraged_2pts("oA0Po-sd_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(-1);
  jvec K_P5P5 =load_chaveraged_2pts("oPPo-ss_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(1);
  jvec Ks_P5P5=load_chaveraged_2pts("oPPo-sd_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(1);
  
  //define the ratios
  jvec ratio_P5P5=Ks_P5P5/K_P5P5;
  jvec ratio_A0P5=Ks_A0P5/K_A0P5;
  
  //fit of P5P5 for K
  jack A_P5P5(njack),SL_P5P5(njack),C_P5P5(njack),M_P5P5(njack);
  jack_fit_mass_and_ratio_P5P5(A_P5P5,SL_P5P5,C_P5P5,M_P5P5,K_P5P5,ratio_P5P5,12,23);
  
  //fit of A0P5 for K
  jack A_A0P5(njack),SL_A0P5(njack),C_A0P5(njack),M_A0P5(njack);
  jack_fit_mass_and_ratio_A0P5(A_A0P5,SL_A0P5,C_A0P5,M_A0P5,K_A0P5,ratio_A0P5,12,23);
  
  //fit of P5P5 of Pi
  jack Mpi=constant_fit(effective_mass(Pi_P5P5),12,23);
  
  //calculate mu-md
  jack dml_P5P5=dM_K_fis/SL_P5P5/Zp390;
  jack dml_A0P5=dM_K_fis/SL_A0P5/Zp390;
  //calculate mu,md
  jack mu_P5P5=ml-dml_P5P5*0.5,md_P5P5=ml+dml_P5P5*0.5;
  jack mu_A0P5=ml-dml_A0P5*0.5,md_A0P5=ml+dml_A0P5*0.5;
  
  //calculate fK
  jack fK_P5P5=(mass[0]+mass[1])*sqrt(C_P5P5)/(M_P5P5*sinh(M_P5P5));
  jack fK_A0P5=C_A0P5*Zv390/sqrt(C_P5P5);
  //cout<<(mass[0]+mass[1])<<" "<<sqrt(C_P5P5)<<" "<<M_P5P5<<" "<<sinh(M_P5P5)<<endl;
  //calculate delta_fK/fK/delta_m
  jack dfK_fr_afK_2dm=A_A0P5-0.5*A_P5P5-0.5*(1/M_P5P5-TH)*SL_A0P5;
  jack dfK_fr_afK_2dm_WI=-1/(mass[0]+mass[1])+0.5*(A_P5P5+(TH-3.0/M_P5P5)*SL_P5P5);
  //physical units
  jack dfK_fr_fK=dfK_fr_afK_2dm*dml_P5P5*Zp390*a390;
  jack dfK_fr_fK_WI=dfK_fr_afK_2dm_WI*dml_A0P5*Zp390*a390;
  
  cout<<"Delta MK2: "<<Mpi*Mpi<<" "<<2*M_P5P5*SL_P5P5<<endl;
  cout<<" md-mu P5P5: "<<dml_P5P5<<", A0P5: "<<dml_A0P5<<endl;
  cout<<" mu/md P5P5: "<<mu_P5P5/md_P5P5<<", A0P5: "<<mu_A0P5/md_A0P5<<endl;
  cout<<endl;
  cout<<"Mass P5P5: "<<M_P5P5<<" = "<<M_P5P5/a390<<" GeV"<<endl;
  cout<<"Slope P5P5: "<<SL_P5P5<<endl;
  cout<<"A P5P5: "<<A_P5P5<<endl;
  cout<<endl;
  cout<<"Mass A0P5: "<<M_A0P5<<" = "<<M_A0P5/a390<<" GeV"<<endl;
  cout<<"Slope A0P5: "<<SL_A0P5<<endl;
  cout<<"A A0P5: "<<A_A0P5<<endl;
  cout<<endl;
  cout<<"fK (def): "<<fK_A0P5<<" = "<<fK_A0P5/a390<<" GeV"<<endl;
  cout<<"fK (WI): "<<fK_P5P5<<" = "<<fK_P5P5/a390<<" GeV"<<endl;
  cout<<endl;
  cout<<"1/a*dfK/fK/2dm (def): "<<dfK_fr_afK_2dm<<endl;
  cout<<"1/a*dfK/fK/2dm (WI): "<<dfK_fr_afK_2dm_WI<<endl;
  cout<<endl;
  cout<<"dfK/fK (def): "<<dfK_fr_fK<<endl;
  cout<<"dfK/fK (WI): "<<dfK_fr_fK_WI<<endl;
  
  (Mpi*Mpi).write_to_binfile("slope_results");
  (2*M_P5P5*SL_P5P5).append_to_binfile("slope_results");
  (2*M_A0P5*SL_A0P5).append_to_binfile("slope_results");
  dfK_fr_fK.append_to_binfile("slope_results");
  dfK_fr_fK_WI.append_to_binfile("slope_results");
  
  return 0;
}
