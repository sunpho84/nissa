#include <include.h>
#include <iostream>

using namespace std;

const int LL=0,LS=1,SL=2,SS=3;
//int tmin[4]={16,6,10,6};
int tmin[4]={18,6,6,12};
jvec c[4];

int nm1=19;
int ijack_fit=0;
int T=48,TH=24;
int njack=10;
double mass[19]={0.006400,0.015000,0.018000,0.022000,0.027000,0.204900,0.230000,0.258200,0.289800,0.325300,0.365100,0.409800,0.460000,0.516300,0.579600,0.650500,0.730200,0.819600,0.919900};

double fun_fit(double Z1,double Z2,double M,int t)
{
  return Z1*Z2*exp(-M*TH)*cosh(M*(TH-t))/M;
}

void ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double ZL0=p[0];
  double M0=p[1];
  double ZL=p[2];
  double ZS=p[3];
  double M=p[4];
  for(int t=0;t<=TH;t++)
    {
      if(t>=tmin[LL] && t<=23) ch+=sqr((c[LL][t].data[ijack_fit]-fun_fit(ZL0,ZL0,M0,t))/c[LL][t].err());
      if(t>=tmin[LL]) ch+=sqr((c[LL][t].data[ijack_fit]-fun_fit(ZL,ZL,M,t))/c[LL][t].err());
      if(t>=tmin[LS]) ch+=sqr((c[LS][t].data[ijack_fit]-fun_fit(ZL,ZS,M,t))/c[LS][t].err());
      if(t>=tmin[SL]) ch+=sqr((c[SL][t].data[ijack_fit]-fun_fit(ZL,ZS,M,t))/c[SL][t].err());
      if(t>=tmin[SS]) ch+=sqr((c[SS][t].data[ijack_fit]-fun_fit(ZS,ZS,M,t))/c[SS][t].err());
    }
}

int icombo(int im1,int im2,int r1,int r2)
{return r1+2*im1+nm1*2*(im2*2+r2);}

int main()
{
  const char tag[4][4]={"LL","LS","SL","SS"};
  
  TMinuit minu;
  minu.SetFCN(ch2);
  minu.SetPrintLevel(-1);
  
  int im2=0;
  jvec a(14,njack),b(14,njack);
  int icorr=0;
  for(int im1=5;im1<19;im1++)
    {
      cout<<"Mass: "<<im1<<" + "<<mass[im1]<<" "<<mass[im2]<<endl;
      cout<<"-----------------------------------"<<endl;
      
      jack M_eff[4];
      
      for(int ism=0;ism<4;ism++)
	{
	  c[ism]=(jvec_load(combine("/tmp/%s",tag[ism]).c_str(),T,njack,icombo(im1,im2,0,0))+jvec_load(combine("/tmp/%s",tag[ism]).c_str(),T,njack,icombo(im1,im2,1,1)))/2;
	  c[ism].print_to_file("/tmp/c");
	  c[ism]=c[ism].simmetrized(1);
	  c[ism].print_to_file(combine("/tmp/csimm_%0.4f_%s",mass[im1],tag[ism]).c_str());
	  jvec e=effective_mass(c[ism]);
	  e.print_to_file(combine("/tmp/eff_%0.04f_%s",mass[im1],tag[ism]).c_str());
	  M_eff[ism]=constant_fit(e,tmin[ism],23);
	  cout<<tag[ism]<<" "<<M_eff[ism]<<endl;
	}
      
      jack ZL0(njack),ZL(njack),ZS(njack),M(njack),M0(njack);
      for(ijack_fit=0;ijack_fit<njack+1;ijack_fit++)
	{
	  double est=M_eff[0][ijack_fit];
	  if(isnan(est)) est=M_eff[2][ijack_fit];
	  
	  minu.DefineParameter(0,"ZL0",0.4,0.001,0,100);
	  minu.DefineParameter(1,"M0",est,0.001,0,0);
	  minu.DefineParameter(2,"ZL",0.4,0.001,0,100);
	  minu.DefineParameter(3,"ZS",0.02,0.001,0,100);
	  minu.DefineParameter(4,"M",M_eff[2][ijack_fit],0.001,0,0);
	  
	  double dum;
	  
	  if(ijack_fit==0) for(int ncall=0;ncall<100;ncall++)
	  minu.mnseek();
	  minu.mnsimp();
	  minu.Migrad();
	  minu.mnmnos();
	  minu.GetParameter(0,ZL0.data[ijack_fit],dum);
	  minu.GetParameter(1,M0.data[ijack_fit],dum);
	  minu.GetParameter(2,ZL.data[ijack_fit],dum);
	  minu.GetParameter(3,ZS.data[ijack_fit],dum);
	  minu.GetParameter(4,M.data[ijack_fit],dum);
	}
      
      cout<<endl;
      cout<<"ZL0 "<<ZL0<<", ZL "<<ZL<<endl; 
      cout<<"M0 "<<M0<<", M "<<M<<endl;
      cout<<"ZS "<<ZS<<endl;
      cout<<"ZL0^2 "<<ZL0*ZL0<<endl;
      
      double lat=1/2.3286;
      double amc=mass[im1];
      double ams=mass[im2];
      cout<<endl;
      cout<<amc/lat/0.44<<" af: "<<ZL*(amc+ams)/(M*sinh(M))<<endl;
      cout<<amc/lat/0.44<<" af0: "<<ZL0*(amc+ams)/(M0*sinh(M0))<<endl;
      cout<<amc/lat/0.44<<" f: "<<ZL*(amc+ams)/(M*sinh(M))/lat<<endl;
      cout<<amc/lat/0.44<<" f0: "<<ZL0*(amc+ams)/(M0*sinh(M0))/lat<<endl;
      
      cout<<"-----------------------------------"<<endl<<endl<<endl;
      
      a[icorr]=ZL*(amc+ams)/(M*sinh(M))/lat;
      b[icorr++]=M/lat;
    }
  
  a.write_to_binfile("fB");
  b.write_to_binfile("mB");
  cout<<(a*sqrt(b));
  
  return 0;
}
