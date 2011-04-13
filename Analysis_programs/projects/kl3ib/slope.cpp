#include "include.h"
#include "kl3_common.cpp"
int njack=10;

double fun_A0P5(double Z_A0P5,double M_A0P5,int t){return Z_A0P5*exp(-M_A0P5*TH)*sinh(M_A0P5*(t-TH));}
double fun_P5P5(double Z_P5P5,double M_P5P5,int t){return Z_P5P5*exp(-M_P5P5*TH)*cosh(M_P5P5*(TH-t))/M_P5P5;}
double fun_ratio_P5P5(double A,double SL,double M,int t){return A+SL*(t-TH)*tanh(M*(t-TH));}
double fun_ratio_A0P5(double A,double SL,double M,int t){return A+SL*(TH-t)/tanh(M*(TH-t));}

/////////////////  fit of slope and ratio of A0P5  /////////////////

double (*fit_fun_K)(double,double,int),(*fit_fun_ratio)(double,double,double,int);
int tmin_fit,tmax_fit;
double *K_fit,*ratio_fit;
double *dK_fit,*dratio_fit;

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
  return (read_two_points(path,T,njack,nmoms,nmass,im1,im2,ik1,ik2,0,0,ri)+
	  read_two_points(path,T,njack,nmoms,nmass,im1,im2,ik1,ik2,1,1,ri))*
    (-0.5/(L*L*L));
}

int main()
{
  read_input();
  
  int im1=0,im2=1; //quello raddoppiato e' il primo!
  int ik1=0,ik2=0;
  int ri=0;
  
  jvec K_A0P5 =load_chaveraged_2pts("oA0Po-ss_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(-1);
  jvec Ks_A0P5=load_chaveraged_2pts("oA0Po-sd_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(-1);
  jvec K_P5P5 =load_chaveraged_2pts("oPPo-ss_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(1);
  jvec Ks_P5P5=load_chaveraged_2pts("oPPo-sd_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(1);
  
  jvec ratio_P5P5=Ks_P5P5/K_P5P5;
  jvec ratio_A0P5=Ks_A0P5/K_A0P5;
  
  jack A_P5P5(njack),SL_P5P5(njack),C_P5P5(njack),M_P5P5(njack);
  jack_fit_mass_and_ratio_P5P5(A_P5P5,SL_P5P5,C_P5P5,M_P5P5,K_P5P5,ratio_P5P5,13,23);

  jack A_A0P5(njack),SL_A0P5(njack),C_A0P5(njack),M_A0P5(njack);
  jack_fit_mass_and_ratio_A0P5(A_A0P5,SL_A0P5,C_A0P5,M_A0P5,K_A0P5,ratio_A0P5,13,23);
  
  cout<<"Mass P5P5: "<<M_P5P5<<endl;
  cout<<"Slope P5P5: "<<SL_P5P5<<endl;
  cout<<"A P5P5: "<<A_P5P5<<endl;
  cout<<endl;
  cout<<"Mass A0P5: "<<M_A0P5<<endl;
  cout<<"Slope A0P5: "<<SL_A0P5<<endl;
  cout<<"A A0P5: "<<A_A0P5<<endl;
  cout<<endl;
  
  jvec rest_ratio_P5P5(ratio_P5P5);
  for(int t=0;t<=TH;t++)
    for(int ijack=0;ijack<njack+1;ijack++)
      rest_ratio_P5P5.data[t].data[ijack]-=fun_ratio_P5P5(A_P5P5[ijack],SL_P5P5[ijack],M_P5P5[ijack],t);
  
  ofstream fuf("/tmp/fuf_jack");
  //cout<<ratio_P5P5.data[0].data<<" "<<rest_ratio_P5P5.data[0].data<<endl;
  fuf<<K_P5P5;
  
  return 0;
}
