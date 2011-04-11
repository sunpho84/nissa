#pragma once

#include <fun.c>

//return effective mass
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
