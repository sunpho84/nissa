#include <iostream>
#include <math.h>

#include "quad.cpp"

using namespace std;

//using gnu
__float128 __float128_sqrt_with_guess(__float128 in,__float128 guess)
{
 return guess+(in-guess*guess)/(2*guess);
}

double conv_double(__float128 in)
{return (double)in;}

__float128 conv___float128(float128 in)
{return (__float128)in.a+(__float128)in.b;}

__float128 __float128_sqrt(__float128 in)
{
  __float128 out;
  __float128 guess=sqrt(in);
  __float128 err;

  do
    {
      out=__float128_sqrt_with_guess(in,guess);
      err=out-guess;
      guess=out;
      //cout<<"errf: "<<(double)err<<endl;
    }
  while(err>1.e-60||err<-1.e-60);
  
  return out;
}

int main()
{
  float128 sqrt3=float128_sqrt(3);
  cout<<"float128 - double:      "<<conv_double(  float128_sqrt(  (float128)3)-sqrt(3))<<endl;
  cout<<"__float128 - double:    "<<conv_double(__float128_sqrt((__float128)3)-sqrt(3))<<endl;
  cout<<"__float128 vs float128: "<<conv_double((float128)__float128_sqrt((__float128)3)-float128_sqrt((float128)3))<<endl;
  cout<<"__float128 vs float128: "<<conv_double(__float128_sqrt((__float128)3)-conv___float128(float128_sqrt((float128)3)))<<endl;

  /*
  __float128 a=float128_sqrt(3),b=1;
  float128 sa=sqrt(float128(3)),sb=1;
  
  for(int i=0;i<800;i++)
    {
      cout<<i<<"\t"<<conv_double(b)<<"\t"<<conv_double((((a+b)-a))-b)<<"\t\t"<<conv_double((((sa+sb)-sa))-sb)<<endl;
      b/=2;
      sb/=2;
    }
  */
  __float128 one_third=1/(__float128)31;
  float128 one_third_sanfo=((float128)1/(float128)31);
  cout<<(double)(one_third-1.0/31)<<endl;
  cout<<(double)(one_third-conv___float128(one_third_sanfo))<<endl;
  return 0;
}
