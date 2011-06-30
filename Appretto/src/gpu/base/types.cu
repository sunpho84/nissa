#pragma once

//single prec types
typedef float complef[2];

typedef complef spinf[4];
typedef complef colorf[3];

typedef colorf su3f[3];
typedef su3f quad_su3f[4];

typedef float su3c[8];


//double prec types
typedef double compled[2];

typedef compled spind[4];
typedef compled colord[3];

typedef colord su3d[3];
typedef su3d quad_su3d[4];


//The sum of two complex number
void complef_summ(complef a,complef b,complef c)
{
  a[0]=b[0]+c[0];
  a[1]=b[1]+c[1];
}
void complef_subt(complef a,complef b,complef c)
{
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
}

//prod with real
void complef_prod_with_real(complef a,complef b,float c)
{
  a[0]=b[0]*c;
  a[1]=b[1]*c;
}

//Summ to the output the product of two complex number
//it is assumed that a!=b and a!=c

void complef_summ_the_prod(complef a,complef b,complef c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
}
void complef_summ_the_conj2_prod(complef a,complef b,complef c)
{
  a[0]+=b[0]*c[0]+b[1]*c[1];
  a[1]+=-b[0]*c[1]+b[1]*c[0];
}
void complef_summ_the_conj1_prod(complef a,complef b,complef c)
{
  complef_summ_the_conj2_prod(a,c,b);
}
void complef_summ_the_conj_conj_prod(complef a,complef b,complef c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=-b[0]*c[1]-b[1]*c[0];
}
void complef_subt_the_prod(complef a,complef b,complef c)
{
  a[0]-=b[0]*c[0]-b[1]*c[1];
  a[1]-=b[0]*c[1]+b[1]*c[0];
}
void complef_subt_the_conj2_prod(complef a,complef b,complef c)
{
  a[0]-=b[0]*c[0]+b[1]*c[1];
  a[1]-=-b[0]*c[1]+b[1]*c[0];
}
void complef_subt_the_conj1_prod(complef a,complef b,complef c)
{
  complef_subt_the_conj2_prod(a,c,b);
}
void complef_subt_the_conj_conj_prod(complef a,complef b,complef c)
{
  a[0]-=b[0]*c[0]-b[1]*c[1];
  a[1]-=-b[0]*c[1]-b[1]*c[0];
}

//The product of two complex number
void unsafe_complef_prod(complef a,complef b,complef c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
}

//The product of a complex number by the conjugate of the second
void unsafe_complef_conj2_prod(complef a,complef b,complef c)
{
  a[0]=b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
}
void unsafe_complef_conj1_prod(complef a,complef b,complef c)
{
  unsafe_complef_conj2_prod(a,c,b);
}

//The product of the conjugate of two complex numbers
void unsafe_complef_conj_conj_prod(complef a,complef b,complef c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=-b[0]*c[1]-b[1]*c[0];
}

//The product of two complex numbers
void safe_complef_prod(complef a,complef b,complef c)
{
  float tmp=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}

//The product of a complex number by the conjugate of the second
void safe_complef_conj2_prod(complef a,complef b,complef c)
{
  float tmp=a[0]=b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}
void safe_complef_conj1_prod(complef a,complef b,complef c)
{
  safe_complef_conj2_prod(a,c,b);
}

//complex prod real
void complef_prod_real(complef a,complef b,float c)
{
  a[0]=b[0]*c;
  a[1]=b[1]*c;
}

//comple prod i
void safe_complef_prod_i(complef a,complef b)
{
  float temp=b[0];
  a[0]=-b[1];
  a[1]=temp;
}
void assign_complef_prod_i(complef a){safe_complef_prod_i(a,a);}

//complex prod -i
void safe_complef_prod_minus_i(complef a,complef b)
{
  float temp=b[0];
  a[0]=b[1];
  a[1]=-temp;
}
void assign_complef_prod_minus_i(complef a){safe_complef_prod_minus_i(a,a);}

//reciprocal of a complex
void complef_reciprocal(complef rec,complef c)
{
  float module=c[0]*c[0]+c[1]*c[1];
  
  rec[0]=c[0]/module;
  rec[1]=-c[1]/module;
}

//power of a complex
void complef_pow(complef res,complef base,float exp)
{
  float module=pow(base[0]*base[0]+base[1]*base[1],exp/2);
  float anomaly=atan2(base[1],base[0])*exp;

  res[0]=module*cos(anomaly);
  res[1]=module*sin(anomaly);
}

//the real amd imaginary unit
complef ONEf={1,0};
complef If={0,1};

//pi
const float PIf=3.14159265358979323846f;
