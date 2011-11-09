#pragma once

//Assign the conj
void complex_conj(complex a,complex b)
{
  a[0]=b[0];
  a[1]=-b[1];
}
//The sum of two complex number
void complex_summ(complex a,complex b,complex c)
{
  a[0]=b[0]+c[0];
  a[1]=b[1]+c[1];
}
void complex_summ_conj2(complex a,complex b,complex c)
{
  a[0]=b[0]+c[0];
  a[1]=b[1]-c[1];
}
void complex_summ_conj1(complex a,complex b,complex c)
{complex_summ_conj2(a,c,b);}
void complex_subt(complex a,complex b,complex c)
{
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
}
void complex_subt_conj2(complex a,complex b,complex c)
{
  a[0]=b[0]-c[0];
  a[1]=b[1]+c[1];
}
void complex_subt_conj1(complex a,complex b,complex c)
{
  a[0]=b[0]-c[0];
  a[1]=-b[1]-c[1];
}
void complex_summassign(complex a,complex b)
{
  a[0]+=b[0];
  a[1]+=b[1];
}
void complex_subtassign(complex a,complex b)
{
  a[0]-=b[0];
  a[1]-=b[1];
}

//prod with real
void complex_prod_with_real(complex a,complex b,double c)
{
  a[0]=b[0]*c;
  a[1]=b[1]*c;
}

//Summ to the output the product of two complex number
//it is assumed that a!=b and a!=c
void complex_summ_the_prod(complex a,complex b,complex c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
}
void complex_subt_the_prod(complex a,complex b,complex c)
{
  a[0]-=b[0]*c[0]-b[1]*c[1];
  a[1]-=b[0]*c[1]+b[1]*c[0];
}
void complex_summ_the_conj2_prod(complex a,complex b,complex c)
{
  a[0]+=+b[0]*c[0]+b[1]*c[1];
  a[1]+=-b[0]*c[1]+b[1]*c[0];
}
void complex_summ_the_conj1_prod(complex a,complex b,complex c)
{complex_summ_the_conj2_prod(a,c,b);}
void complex_summ_the_conj_conj_prod(complex a,complex b,complex c)
{
  a[0]+=+b[0]*c[0]-b[1]*c[1];
  a[1]+=-b[0]*c[1]-b[1]*c[0];
}
void complex_subt_the_conj2_prod(complex a,complex b,complex c)
{
  a[0]-=+b[0]*c[0]+b[1]*c[1];
  a[1]-=-b[0]*c[1]+b[1]*c[0];
}
void complex_subt_the_conj1_prod(complex a,complex b,complex c)
{complex_subt_the_conj2_prod(a,c,b);}
void complex_subt_the_conj_conj_prod(complex a,complex b,complex c)
{
  a[0]-=+b[0]*c[0]-b[1]*c[1];
  a[1]-=-b[0]*c[1]-b[1]*c[0];
}

//The product of two complex number
void unsafe_complex_prod(complex a,complex b,complex c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
}

//The product of a complex number by the conjugate of the second
void unsafe_complex_conj2_prod(complex a,complex b,complex c)
{
  a[0]=+b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
}
void unsafe_complex_conj1_prod(complex a,complex b,complex c)
{unsafe_complex_conj2_prod(a,c,b);}

//The product of the conjugate of two complex numbers
void unsafe_complex_conj_conj_prod(complex a,complex b,complex c)
{
  a[0]=+b[0]*c[0]-b[1]*c[1];
  a[1]=-b[0]*c[1]-b[1]*c[0];
}

//The product of two complex number
void safe_complex_prod(complex a,complex b,complex c)
{
  double tmp=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}

//The product of a complex number by the conjugate of the second
void safe_complex_conj2_prod(complex a,complex b,complex c)
{
  double tmp=b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}
void safe_complex_conj1_prod(complex a,complex b,complex c)
{safe_complex_conj2_prod(a,c,b);}

//complex prod real
void complex_prod_real(complex a,complex b,double c)
{
  a[0]=b[0]*c;
  a[1]=b[1]*c;
}

//complex prod i
void safe_complex_prod_i(complex a,complex b)
{
  double tmp=b[0];
  a[0]=-b[1];
  a[1]=tmp;
}
void assign_complex_prod_i(complex a)
{safe_complex_prod_i(a,a);}

//complex prod -i
void safe_complex_prod_minus_i(complex a,complex b)
{
  double tmp=b[0];
  a[0]=b[1];
  a[1]=-tmp;
}
void assign_complex_prod_minus_i(complex a)
{safe_complex_prod_minus_i(a,a);}
void complex_summ_the_prod_i(complex a,complex b,complex c)
{
  a[1]+=b[0]*c[0]-b[1]*c[1];
  a[0]-=b[0]*c[1]+b[1]*c[0];
}
void complex_subt_the_prod_i(complex a,complex b,complex c)
{
  a[1]-=b[0]*c[0]-b[1]*c[1];
  a[0]+=b[0]*c[1]+b[1]*c[0];
}
void complex_summ_the_conj2_prod_i(complex a,complex b,complex c)
{
  a[1]+=+b[0]*c[0]+b[1]*c[1];
  a[0]-=-b[0]*c[1]+b[1]*c[0];
}
void complex_summ_the_conj1_prod_i(complex a,complex b,complex c)
{complex_summ_the_conj2_prod(a,c,b);}
void complex_summ_the_conj_conj_prod_i(complex a,complex b,complex c)
{
  a[1]+=+b[0]*c[0]-b[1]*c[1];
  a[0]-=-b[0]*c[1]-b[1]*c[0];
}
void complex_subt_the_conj2_prod_i(complex a,complex b,complex c)
{
  a[1]-=+b[0]*c[0]+b[1]*c[1];
  a[0]+=-b[0]*c[1]+b[1]*c[0];
}
void complex_subt_the_conj1_prod_i(complex a,complex b,complex c)
{complex_subt_the_conj2_prod(a,c,b);}
void complex_subt_the_conj_conj_prod_i(complex a,complex b,complex c)
{
  a[1]-=+b[0]*c[0]-b[1]*c[1];
  a[0]+=-b[0]*c[1]-b[1]*c[0];
}
//squared norm
double squared_complex_norm(complex c)
{return c[0]*c[0]+c[1]*c[1];}

//reciprocal of a complex
void complex_reciprocal(complex rec,complex c)
{
  double module=c[0]*c[0]+c[1]*c[1];
  
  rec[0]=c[0]/module;
  rec[1]=-c[1]/module;
}

//squared root of a complex
void complex_sqrt(complex res,complex base)
{
  double module=sqrt(base[0]*base[0]+base[1]*base[1]);
  double cost=base[0]/module;
  double sinth=sqrt(0.5*(1-cost));
  double costh=sqrt(0.5*(1+cost));
  module=sqrt(module);
  if(base[1]>=0) res[0]=+module*costh;
  else           res[0]=-module*costh;
  res[1]=module*sinth;
}

//power of a complex
void complex_pow(complex res,complex base,double exp)
{
  double module=pow(base[0]*base[0]+base[1]*base[1],exp/2);
  double anomaly=atan2(base[1],base[0])*exp;

  res[0]=module*cos(anomaly);
  res[1]=module*sin(anomaly);
}
