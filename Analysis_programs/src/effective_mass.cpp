#pragma once

//return effective mass

double effective_mass(double y,double yp,int t,int TH)
{
  double rapp=yp/y;
  double meff=-log(rapp);

  double res;
  int iter=0;
  do
    {
      double rapr=cosh(meff*(TH-(t+1)))/cosh(meff*(TH-t));
      double scale=rapr/rapp;
      res=fabs(1-scale)/meff;
      meff*=scale;
      iter++;
    }
  while(res>1.e-8 && iter<1000);
  
  return meff;
}
 //return effective mass
 /*
double fun_cosh(double C,double M,int t,int TH)
{return C*exp(-M*TH)*cosh(M*(TH-t))/M;}
     
double effective_mass(double y,double yp,int t,int TH)
{
  double rapp=yp/y;
  double meff=-log(rapp);

  double res;
  int iter=0;
  do
    {
      double rapr=fun_cosh(1,meff,t+1,TH)/fun_cosh(1,meff,t,TH);
      double scale=rapr/rapp;
      res=fabs(1-scale)/meff;
      meff*=scale;
      iter++;
    }
  while(res>1.e-8 && iter<1000);

  return meff;
}
*/

jack effective_mass(const jack y,const jack yp,int t,int TH)
{
  int njack=y.njack;
  jack m(njack);
  
  for(int ij=0;ij<njack+1;ij++) m.data[ij]=effective_mass(y.data[ij],yp.data[ij],t,TH);
  
  return m;
}

jvec effective_mass(const jvec a)
{
  int TH=a.nel-1;
  int njack=a.njack;
  
  jvec b(TH,njack);
  
  for(int t=0;t<TH;t++) b.data[t]=effective_mass(a.data[t],a.data[t+1],t,TH);
  
  return b;
}
