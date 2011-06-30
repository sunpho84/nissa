#pragma once

//single prec types
typedef float complef[2];

typedef complef spinf[4];
typedef complef colorf[3];

typedef colorf su3f[3];
typedef su3f quad_su3f[4];

typedef float su3c[8];
typedef su3c quad_su3c[4];


//double prec types
typedef double compled[2];

typedef compled spind[4];
typedef compled colord[3];

typedef colord spincolord[4];

typedef colord su3d[3];
typedef su3d quad_su3d[4];

//The sum of two complex number
#define complef_summ(a,b,c) a[0]=b[0]+c[0];a[1]=b[1]+c[1];
#define complef_subt(a,b,c) a[0]=b[0]-c[0];a[1]=b[1]-c[1];
//prod with real
#define complef_prod_with_real(a,b,c) a[0]=b[0]*c;a[1]=b[1]*c;
//Summ to the output the product of two complex number. It is assumed that a!=b and a!=c
#define complef_summ_the_prod(a,b,c) a[0]+=b[0]*c[0]-b[1]*c[1];a[1]+=b[0]*c[1]+b[1]*c[0];
#define complef_summ_the_conj2_prod(a,b,c) a[0]+=b[0]*c[0]+b[1]*c[1];a[1]+=-b[0]*c[1]+b[1]*c[0];
#define complef_summ_the_conj1_prod(a,b,c) complef_summ_the_conj2_prod(a,c,b);
#define complef_summ_the_conj_conj_prod(a,b,c) a[0]+=b[0]*c[0]-b[1]*c[1]; a[1]-=b[0]*c[1]+b[1]*c[0];
#define complef_subt_the_prod(a,b,c) a[0]-=b[0]*c[0]-b[1]*c[1];a[1]-=b[0]*c[1]+b[1]*c[0];
#define complef_subt_the_conj2_prod(a,b,c) a[0]-=b[0]*c[0]+b[1]*c[1];a[1]-=-b[0]*c[1]+b[1]*c[0];
#define complef_subt_the_conj1_prod(a,b,c) complef_subt_the_conj2_prod(a,c,b);
#define complef_subt_the_conj_conj_prod(a,b,c) a[0]-=b[0]*c[0]-b[1]*c[1];a[1]+=b[0]*c[1]+b[1]*c[0];
//The product of two complex number
#define unsafe_complef_prod(a,b,c) a[0]=b[0]*c[0]-b[1]*c[1];a[1]=b[0]*c[1]+b[1]*c[0];
//The product of a complex number by the conjugate of the second
#define unsafe_complef_conj2_prod(a,b,c) a[0]=b[0]*c[0]+b[1]*c[1];a[1]=-b[0]*c[1]+b[1]*c[0];
#define unsafe_complef_conj1_prod(a,b,c) unsafe_complef_conj2_prod(a,c,b);
//The product of the conjugate of two complex numbers
#define unsafe_complef_conj_conj_prod(a,b,c) a[0]=b[0]*c[0]-b[1]*c[1];a[1]=-b[0]*c[1]-b[1]*c[0];
//The product of two complex numbers
#define safe_complef_prod(a,b,c,tmp) tmp=b[0]*c[0]-b[1]*c[1];a[1]=b[0]*c[1]+b[1]*c[0];a[0]=tmp;
//The product of a complex number by the conjugate of the second
#define safe_complef_conj2_prod(a,b,c,tmp) tmp=b[0]*c[0]+b[1]*c[1];a[1]=-b[0]*c[1]+b[1]*c[0];a[0]=tmp;
#define safe_complef_conj1_prod(a,b,c,tmp)safe_complef_conj2_prod(a,c,b,tmp);
#define safe_complef_conj_conj_prod(a,b,c,tmp) tmp=b[0]*c[0]-b[1]*c[1];a[1]=-b[0]*c[1]-b[1]*c[0];a[0]=tmp;
//complex prod real
#define complef_prod_real(a,b,c) a[0]=b[0]*c;a[1]=b[1]*c;
//comple prod i
#define safe_complef_prod_i(a,b,tmp) tmp=b[0];a[0]=-b[1];a[1]=tmp;
#define safe_assign_complef_prod_i(a,tmp) safe_complef_prod_i(a,a,tmp);
//complex prod -i
#define safe_complef_prod_minus_i(a,b,tmp) tmp=b[0];a[0]=b[1];a[1]=-tmp;
#define safe_assign_complef_prod_minus_i(a,tmp) safe_complef_prod_minus_i(a,a,tmp);

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
