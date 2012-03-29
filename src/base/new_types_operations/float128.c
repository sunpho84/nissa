#pragma once

void float_128_uminus(float_128 b,float_128 a)
{
  b[0]=-a[0];
  b[1]=-a[1];
}

void float_128_summ(float_128 c,float_128 a,float_128 b)
{
  double t1=a[0]+b[0];
  double e=t1-a[0];
  double t2=((b[0]-e)+(a[0]-(t1-e)))+a[1]+b[1];
  
  c[0]=t1+t2;
  c[1]=t2-(c[0]-t1);
}
void float_128_summassign(float_128 b,float_128 a)
{float_128_summ(b,b,a);}
void float_128_summ_double(float_128 c,float_128 a,double b)
{
  float_128 d={b,0};
  float_128_summ(c,a,d);
}

void float_128_summassign_double(float_128 a,double b)
{float_128_summ_double(a,a,b);}

void float_128_subt_from_double(float_128 c,double a,float_128 b)
{
  float_128 d;
  float_128_uminus(d,b);
  float_128_summ_double(c,d,a);
}

void float_128_prod(float_128 c,float_128 a,float_128 b)
{
  double split=134217729;
  
  double cona=a[0]*split;
  double conb=b[0]*split;
  
  double a1=cona-(cona-a[0]);
  double b1=conb-(conb-b[0]);
  double a2=a[0]-a1;
  double b2=b[0]-b1;
  
  double c11=a[0]*b[0];
  double c21=a2*b2+(a2*b1+(a1*b2+(a1*b1-c11)));
  
  double c2=a[0]*b[1]+a[1]*b[0];
  
  double t1=c11+c2;
  double e=t1-c11;
  double t2=a[1]*b[1]+((c2-e)+(c11-(t1-e)))+c21;
  
  c[0]=t1+t2;
  c[1]=t2-(c[0]-t1);
}

void float_128_summ_the_prod(float_128 c,float_128 a,float_128 b)
{
  float_128 d;
  float_128_prod(d,a,b);
  float_128_summassign(c,d);
}
