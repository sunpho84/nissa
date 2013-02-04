#include <stdint.h>
#include <string.h>
#include <math.h>

#include "../base/debug.h"
#include "../base/macros.h"
#include "../base/routines.h"
#include "float128.h"

#include "float256bis.cpp"

void float_256_print(float_256 a)
{
  master_printf("(");
  float_128_print((double*)a);
  master_printf(",");
  float_128_print((double*)a+2);
  master_printf(")");
}

void float_256_copy(float_256 b,float_256 a)
{for(int i=0;i<4;i++) b[i]=a[i];}

void float_256_swap(float_256 b,float_256 a)
{
  float_256 c;
  float_256_copy(c,b);
  float_256_copy(b,a);
  float_256_copy(a,c);
}

void float_256_uminus(float_256 b,float_256 a)
{for(int i=0;i<4;i++) b[i]=-a[i];}

void float_256_from_double(float_256 b,double a)
{
  b[0]=a;
  b[1]=b[2]=b[3]=0;
}

void float_256_put_to_zero(float_256 a)
{float_256_from_double(a,0);}

void float_256_summassign(float_256 b,float_256 a)
{float_256_summ(b,b,a);}
void float_256_subt(float_256 c,float_256 a,float_256 b)
{
  float_256 d;
  float_256_uminus(d,b);
  float_256_summ(c,d,a);
}
void float_256_subtassign(float_256 b,float_256 a)
{float_256_subt(b,b,a);}

void float_256_summ_the_prod(float_256 c,float_256 a,float_256 b)
{
  float_256 d;
  float_256_prod(d,a,b);
  float_256_summassign(c,d);
}
void float_256_subt_the_prod(float_256 c,float_256 a,float_256 b)
{
  float_256 d;
  float_256_prod(d,a,b);
  float_256_subtassign(c,d);
}

void float_256_prodassign(float_256 out,float_256 in)
{float_256_prod(out,out,in);}

void float_256_div(float_256 c,float_256 a,float_256 b)
{
  float_256_unr q;

  float_256 r;  
  float_256_copy(r,a);
  
  for(int i=0;i<4;i++)
    {
      q[i]=r[0]/b[0];
      //printf("remaining %d, ",i);
      //float_256_print(r);
      //printf("\n");
      float_256 s;
      float_256_prod_double(s,b,-q[i]);
      float_256_summassign(r,s);
      //printf("summing ");
      //float_256_print(s);
      //printf("\n");
    }
  
  q[4]=r[0]/b[0];
  
  renormalize(c,q);
}

//integer power
void float_256_pow_int(float_256 out,float_256 in,int d)
{
  //printf("taking %d power of ",d);
  //float_256_print(in);
  //printf("\n");
  
  //negative or null case
  if(d<=0)
    {
      float_256_from_double(out,1);
      //if negative
      if(d<0)
	{
	  //compute inv and assign to out
	  float_256 inv;
	  float_256_div(inv,out,in);
	  float_256_copy(out,inv);
	  //multiply the remaining numer of times
	  for(int i=2;i<=-d;i++) float_256_prodassign(out,inv);
	}
    }
  //positive case
  else
    {
      float_256_copy(out,in);
      for(int i=2;i<=d;i++)
	{
	  //printf("iter %d power, ",i);
	  //float_256_print(out);
	  //printf("\n");
	  float_256_prodassign(out,in);
	}
    }
  
  //printf("result, ");
  //float_256_print(out);
  //printf("\n");
}

//frac power
void float_256_pow_int_frac(float_256 out,float_256 ext_in,int n,int d)
{
  //take a copy
  float_256 in;
  float_256_copy(in,ext_in);
  
  //compute by solving out^d=in^n=ref
  float_256 ref;
  float_256_pow_int(ref,in,n);
  
  //let's start from a reasonable approx
  float_256_from_double(out,pow(in[0],(double)n/d));
  //another possible approach
  //float_256_from_128(out,1+(in[0]-1)*r);

  //comput the reciprocal of d
  float_256 recd,tempd;
  float_256_from_double(recd,1);
  float_256_from_double(tempd,d);
  float_256_div(recd,recd,tempd);

  //printf("  recd: ");
  //float_256_print(recd);
  //printf("\n");
  
  //(out+err)^d=in^n -> err=out*rel_err, rel_err=(ref/out^d-1)/d

  float_256 rel_err;
  do
    {
      //printf("  ref: ");
      //float_256_print(ref);
      //printf("\n");
      
      //compute out^d
      float_256 outd;
      float_256_pow_int(outd,out,d);
      
      //printf("  outd: ");
      //float_256_print(outd);
      //printf("\n");
      
      //divide ref by out^d and subtract 1
      float_256_subt(rel_err,ref,outd);

      //printf("  diff: ");
      //float_256_print(rel_err);
      //printf("\n");
      
      float_256_div(rel_err,rel_err,ref);
      
      //printf("rel err: ");
      //float_256_print(rel_err);
      //printf("\n");
      
      //divide by d
      float_256_prod(rel_err,rel_err,recd);
      
      //total err
      float_256 err;
      float_256_prod(err,rel_err,out);
      
      //printf("err: ");
      //float_256_print(err);
      //printf("\n");
      
      float_256_summassign(out,err);
    }
  while(fabs(rel_err[0])>3.e-57);

  //printf("    out: ");
  //float_256_print(out);
  //printf("\n");    
}

//a>b?
int float_256_is_greater(float_256 a,float_256 b)
{return (a[0]>b[0]||(a[0]==b[0]&&(a[1]>b[1]||(a[1]==b[1]&&(a[2]>b[2]||(a[2]==b[2]&&a[3]>b[3]))))));}
int float_256_is_greater(float_256 a,double b)
{return (a[0]>b||(a[0]==b&&(a[1]>0)));}

//a<b?
int float_256_is_smaller(float_256 a,float_256 b)
{return (a[0]<b[0]||(a[0]==b[0]&&(a[1]<b[1]||(a[1]==b[1]&&(a[2]<b[2]||(a[2]==b[2]&&a[3]<b[3]))))));}
int float_256_is_smaller(float_256 a,double b)
{return (a[0]<b||(a[0]==b&&(a[1]<0)));}

void float_256_abs(float_256 a,float_256 b)
{
  float_256 zero;
  float_256_put_to_zero(zero);
  if(float_256_is_greater(b,zero)) float_256_copy(a,b);
  else float_256_uminus(a,b);
}
