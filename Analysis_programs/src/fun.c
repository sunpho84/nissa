#pragma once

double fun_cosh(double C,double M,int t,int TH)
{
  return C*exp(-M*TH)*cosh(M*(TH-t))/M;
}
