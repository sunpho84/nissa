#pragma once

class ran_gen
{
private:
  static const int ran2_ntab=32;
  int ran2_idum,ran2_idum2,ran2_iy;
  int ran2_iv[ran2_ntab];
public:
  ran_gen(int seed);
  int get_int(int max);
  double get_double();
  double get_double(double max);
  double get_gauss(double med,double sigma);
};

ran_gen::ran_gen(int seed)
{
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  
  ran2_idum=seed;
  
  ran2_idum=max(ran2_idum+1,1);
  ran2_idum2=ran2_idum;
  for(int j=ran2_ntab+7;j>=0;j--)
    {
      int k=ran2_idum/iq1;
      ran2_idum=ia1*(ran2_idum-k*iq1)-k*ir1;
      if(ran2_idum<0) ran2_idum+=im1;
      if(j<ran2_ntab) ran2_iv[j]=ran2_idum;
    }
  ran2_iy=ran2_iv[0];
}

//standard ran2 from numerical recipes
double ran_gen::get_double()
{
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/ran2_ntab;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
  double out;
  
  k=ran2_idum/iq1;
  ran2_idum=ia1*(ran2_idum-k*iq1)-k*ir1;
  if(ran2_idum<0) ran2_idum+=im1;
  
  k=ran2_idum2/iq2;
  ran2_idum2=ia2*(ran2_idum2-k*iq2)-k*ir2;
  if(ran2_idum2<0) ran2_idum2+=im2;
      
  j=ran2_iy/ndiv;
  ran2_iy=ran2_iv[j]-ran2_idum2;
  ran2_iv[j]=ran2_idum;
  if(ran2_iy<0) ran2_iy+=imm1;
  
  out=min(am*ran2_iy,rnmx);
  
  return out;
}

double ran_gen::get_double(double max)
{return get_double()*max;}

int ran_gen::get_int(int max)
{return (int)(get_double(max));}

double ran_gen::get_gauss(double med,double sigma)
{
  double q,r,x;
  static double y;
  static bool flag=true;

  if(flag)
    {
      r=sqrt(-2*log(1-get_double())); //This is correct, otherwise may be inf
      q=2*M_PI*get_double();
      
      x=r*cos(q);
      y=r*sin(q);
    }
  else x=y;

  flag=!flag;

  return med+sigma*x;
}

