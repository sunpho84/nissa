#pragma once

typedef double quadf[2];

void quadf_add(quadf c,quadf a,quadf b)
{
  double t1=a[0]+b[0];
  double e=t1-a[0];
  double t2=((b[0]-e)+(a[0]-(t1-e)))+a[1]+b[1];
  
  c[0]=t1+t2;
  c[1]=t2-(c[0]-t1);
}

void quadf_mul(quadf c,quadf a,quadf b)
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


