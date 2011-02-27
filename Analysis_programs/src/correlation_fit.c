#pragma once

#include <TMinuit.h>

jack_vec *_data_fit_mass;
double *_err_data_fit_mass;
int _ijack_fit_mass;
int _tmin_fit_mass,_tmax_fit_mass;

double chi2_fit_mass(double C,double M,int ij)
{
  int T=_data_fit_mass->nel,TH=T/2;
  double ch2=0;
  
  for(int t=0;t<T;t++)
    if((t>=_tmin_fit_mass && t<=_tmax_fit_mass)||(t>=T-_tmax_fit_mass && t<=T-_tmin_fit_mass))
      {
	double add=pow((_data_fit_mass->data[t][_ijack_fit_mass]-fun_cosh(C,M,t,TH))/_err_data_fit_mass[t],2);
	ch2+=add;
      }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void ch2_fit_mass_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double C=p[0];
  double M=p[1];
  
  ch=chi2_fit_mass(C,M,_ijack_fit_mass);
}

//perform the simultaneous fit of the mass
void jack_fit_mass(jack Z,jack M,jack_vec *a,int tmin,int tmax)
{
  _tmin_fit_mass=tmin;
  _tmax_fit_mass=tmax;
  _err_data_fit_mass=(double*)malloc(sizeof(double)*(a->nel));
  int T=a->nel;
  for(int t=0;t<T;t++)
    if((t>=tmin && t<=tmax)||(t>=T-tmax && t<=T-tmin))
      _err_data_fit_mass[t]=jack_error(a->data[t]);
  _data_fit_mass=a;
 
  for(_ijack_fit_mass=0;_ijack_fit_mass<njack+1;_ijack_fit_mass++)
    {
      TMinuit minu;
      
      minu.SetPrintLevel(-1);
      
      int T=a->nel,TH=T/2;
      double meff=-log(a->data[tmin+1][_ijack_fit_mass]/a->data[tmin][_ijack_fit_mass]);
      double cestim=fun_cosh(1,meff,tmax,TH)/a->data[tmax][_ijack_fit_mass];
      
      minu.DefineParameter(0,"Z",cestim,0.0001,0,0);
      minu.DefineParameter(1,"M",meff,0.0001,0,0);
      
      minu.SetFCN(ch2_fit_mass_wr);
      
      minu.Migrad();
      
      double dum;
      minu.GetParameter(0,Z[_ijack_fit_mass],dum);
      minu.GetParameter(1,M[_ijack_fit_mass],dum);
    }

  free(_err_data_fit_mass);
}
