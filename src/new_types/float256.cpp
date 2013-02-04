#include <stdint.h>
#include <string.h>
#include <math.h>

#include "../base/debug.h"
#include "../base/macros.h"
#include "../base/routines.h"
#include "float128.h"

void float_256_print(float_256 a)
{
  master_printf("(");
  float_128_print(a[0]);
  master_printf(",");
  float_128_print(a[1]);
  master_printf(")");
}

double float_256_most_sign_part(float_256 b)
{return b[0][0];}

void float_256_copy(float_256 b,float_256 a)
{
  float_128_copy(b[0],a[0]);
  float_128_copy(b[1],a[1]);
}

void float_256_swap(float_256 b,float_256 a)
{
  float_256 c;
  float_256_copy(c,b);
  float_256_copy(b,a);
  float_256_copy(a,c);
}

void float_256_uminus(float_256 b,float_256 a)
{
  float_128_uminus(b[0],a[0]);
  float_128_uminus(b[1],a[1]);
}

void float_256_from_128(float_256 b,float_128 a)
{
  float_128_copy(b[0],a);
  float_128_put_to_zero(b[1]);
}

void float_256_from_double(float_256 b,double a)
{
  float_128_from_64(b[0],a);
  float_128_put_to_zero(b[1]);
}

void float_256_put_to_zero(float_256 a)
{
  float_128_put_to_zero(a[0]);
  float_128_put_to_zero(a[1]);
}

void float_128_from_256(float_128 a,float_256 b)
{float_128_summ(a,b[0],b[1]);}

void float_256_summ(float_256 c,float_256 a,float_256 b)
{
  printf("summing: ");
  float_256_print(a);
  printf("+");
  float_256_print(b);
  
  float_128 t1,e,t2,r1,r2,r3,r4;
  float_128_summ(t1,a[0],b[0]);
  float_128_subt(e,t1,a[0]);
  float_128_subt(r1,b[0],e);
  float_128_subt(r2,t1,e);
  float_128_subt(r3,a[0],r2);
  float_128_summ(t2,r1,r3);
  float_128_summassign(t2,a[1]);
  float_128_summassign(t2,b[1]);
  float_128_summ(c[0],t1,t2);
  float_128_subt(r4,c[0],t1);
  float_128_subt(c[1],t2,r4);

  printf("\n = ");
  float_256_print(c);
  printf("\n");
}

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

void float_256_summ_128(float_256 c,float_256 a,float_128 b)
{
  float_256 temp;
  float_256_from_128(temp,b);
  float_256_summ(c,a,temp);
}
void float_256_128_summ_128(float_256 c,float_128 a,float_128 b)
{
  float_256 temp1;
  float_256_from_128(temp1,a);
  float_256 temp2;
  float_256_from_128(temp2,b);
  float_256_summ(c,temp1,temp2);
}
void float_256_summassign_128(float_256 b,float_128 a)
{float_256_summ_128(b,b,a);}
void float_256_subt_from_128(float_256 c,float_128 a,float_256 b)
{
  float_256 d;
  float_256_uminus(d,b);
  float_256_summ_128(c,d,a);
}

void float_128_equi(double &a1,double &a2,float_128 a)
{
  double split=134217729;

  double cona=a[0]*split;

  double t1=cona-(cona-a[0]);
  double t2=a[0]-t1;
  
  a1=t1;
  a2=t2;
}

void float_256_prod(float_256 c,float_256 a,float_256 b)
{
  float_128 split={18014398509481984,1};

  float_128 cona,conb;
  float_128_prod(cona,a[0],split);
  float_128_prod(conb,b[0],split);

  float_128 aa,bb;
  float_128_subt(aa,cona,a[0]);
  float_128_subt(bb,conb,b[0]);
  float_128 a1,b1;
  float_128_subt(a1,cona,aa);
  float_128_subt(b1,conb,bb);
  float_128 a2,b2;
  float_128_subt(a2,a[0],a1);
  float_128_subt(b2,b[0],b1);
  
  float_128_equi(a1[0],a1[1],a1);
  float_128_equi(a2[0],a2[1],a2);
  float_128_equi(b1[0],b1[1],b1);
  float_128_equi(b2[0],b2[1],b2);
  
  float_128 c11;
  float_128_prod(c11,a[0],b[0]);
  float_128 c21,a1b1,t,a1b2,a1b2t,a2b1,a2b1a1b2t;
  float_128_prod(c21,a2,b2);
  float_128_prod(a1b1,a1,b1);
  float_128_subt(t,a1b1,c11);
  float_128_prod(a1b2,a1,b2);
  float_128_summ(a1b2t,a1b2,t);
  float_128_prod(a2b1,a2,b1);
  float_128_summ(a2b1a1b2t,a2b1,a1b2t);
  float_128_summassign(c21,a2b1a1b2t);

  float_128 c2;
  float_128_prod(c2,a[0],b[1]);
  float_128_summ_the_prod(c2,a[1],b[0]);

  float_128 t1,e,t2;
  float_128_summ(t1,c11,c2);
  float_128_subt(e,t1,c11);
  float_128_prod(t2,a[1],b[1]);
  float_128 c2e,t1e,c11t1e,c2ec11t1e;
  float_128_subt(c2e,c2,e);
  float_128_subt(t1e,t1,e);
  float_128_subt(c11t1e,c11,t1e);
  float_128_summ(c2ec11t1e,c2e,c11t1e);
  float_128_summassign(t2,c2ec11t1e);
  float_128_summassign(t2,c21);

  float_128_summ(c[0],t1,t2);
  float_128 c0mt1;
  float_128_subt(c0mt1,c[0],t1);
  float_128_subt(c[1],t2,c0mt1);
}

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

void float_256_summ_the_128_prod(float_256 c,float_128 a,float_128 b)
{
  float_256 temp1;
  float_256_from_128(temp1,a);
  float_256 temp2;
  float_256_from_128(temp2,b);
  float_256_summ_the_prod(c,temp1,temp2);
}

void float_256_prod_128(float_256 c,float_256 a,float_128 b)
{
  float_256 temp;
  float_256_from_128(temp,b);
  float_256_prod(c,a,temp);
}

void float_128_prod_256(float_256 c,float_128 a,float_256 b)
{float_256_prod_128(c,b,a);}
void float_summ_the_128_prod_256(float_256 c,float_128 a,float_256 b)
{
  float_256 d;
  float_128_prod_256(d,a,b);
  float_256_summassign(c,d);
}
void float_subt_the_128_prod_256(float_256 c,float_128 a,float_256 b)
{
  float_256 d;
  float_128_prod_256(d,a,b);
  float_256_subtassign(c,d);
}
void float_256_prodassign(float_256 out,float_256 in)
{float_256_prod(out,out,in);}

//divide two float_256
void float_256_div(float_256 div,float_256 ext_a,float_256 ext_b)
{
  float_256 a,b;
  float_256_copy(a,ext_a);
  float_256_copy(b,ext_b);

  //compute dividing factor
  float_128 c,d;
  float_128_from_64(c,1);
  float_128_summ(d,b[0],b[1]);
  float_128_div(c,c,d);
  //compute approx div
  float_256 div1;
  float_256_prod_128(div1,a,c);
  //compute first remaining
  float_256 rem;
  float_256_prod(rem,div1,b);
  float_256_subt(rem,a,rem);
  //compute the second piece
  float_256 div2;
  float_256_prod_128(div2,rem,c);
  //summ the two pieces
  float_256 div12;
  float_256_summ(div12,div1,div2);
  //second remaining
  float_256_prod(rem,div12,b);
  float_256_subt(rem,a,rem);
  //compute the third piece
  float_256 div3;
  float_256_prod_128(div3,rem,c);
  //summ the two pieces
  float_256_summ(div,div12,div3);
  //third remaining
  float_256_prod(rem,div,b);
  float_256_subt(rem,a,rem);
  printf("third remaining: ");
  float_256_print(rem);
  printf("\n");
}

//divide a float_256 by a double
void float_256_div_128(float_256 div,float_256 a,float_128 b)
{
  float_256 temp;
  float_256_from_128(temp,b);
  float_256_div(div,a,temp);
}

//integer power
void float_256_pow_int(float_256 out,float_256 in,int d)
{
  printf("taking %d power of ",d);
  float_256_print(in);
  printf("\n");
  
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
	  printf("iter %d power, ",i);
	  float_256_print(out);
	  printf("\n");
	  float_256_prodassign(out,in);
	}
    }
  
  printf("result, ");
  float_256_print(out);
  printf("\n");
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
  float_128_pow_int_frac(out[0],in[0],n,d);
  float_128_put_to_zero(out[1]);
  //another possible approach
  //float_256_from_128(out,1+(in[0]-1)*r);

  //comput the reciprocal of d
  float_256 recd,tempd;
  float_256_from_double(recd,1);
  float_256_from_double(tempd,d);
  float_256_div(recd,recd,tempd);
  
  //(out+err)^d=in^n -> err=out*rel_err, rel_err=(ref/out^d-1)/d

  float_256 rel_err;
  do
    {
      //compute out^d
      printf(" ref: ");
      float_256_print(ref);
      printf("\n");
      
      float_256 outd;
      float_256_pow_int(outd,out,d);
      
      //divide ref by out^d and subtract 1
      float_256_div(rel_err,ref,outd);
      float_256 mone;
      float_256_from_double(mone,-1);
      float_256_summassign(rel_err,mone);
      //divide by d
      float_256_prod(rel_err,rel_err,recd);
      
      printf("outd: ");
      float_256_print(outd);
      printf("\n");
      
      printf("rel err: ");
      float_256_print(rel_err);
      printf("\n");
      
      //total err
      float_256 err;
      float_256_prod(err,rel_err,out);
      
      printf("err: ");
      float_256_print(err);
      printf("\n");
      
      float_256_summassign(out,err);
    }
  while(fabs(float_256_most_sign_part(rel_err))>3.e-57);
}

//a>b?
int float_256_is_greater(float_256 a,float_256 b)
{
  if(float_128_is_greater(a[0],b[0])) return true;
  if(float_128_is_greater(b[0],a[0])) return false;
  
  if(float_128_is_greater(a[1],b[1])) return true;
		  
  return false;
}

//a<b?
int float_256_is_smaller(float_256 a,float_256 b)
{return float_256_is_greater(b,a);}

void float_256_abs(float_256 a,float_256 b)
{
  float_256 zero;
  float_256_put_to_zero(zero);
  if(float_256_is_greater(b,zero)) float_256_copy(a,b);
  else float_256_uminus(a,b);
}
