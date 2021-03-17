#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "float_128.hpp"
#include "float_256.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  void float_256_print(float_256 a)
  {
    crash("");
    // printf("(");
    // float_128_print((double*)a);
    // printf(",");
    // float_128_print((double*)a+2);
    // printf(")");
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
  
  void renormalize(float_256,double,double,double,double);
  void float_256_from_double(float_256 b,double a)
  {renormalize(b,a,0,0,0);}
  
  void float_256_put_to_zero(float_256 a)
  {float_256_from_double(a,0);}
  
  //a==b?
  int float_256_is_equal(float_256 a,float_256 b)
  {return a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2]&&a[3]==b[3];}
  
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
  
  void quick_two_summ(double &s,double &e,double a,double b)
  {
    s=a+b;
    e=b-(s-a);
  }
  
  void two_summ(double &s,double &e,double a,double b)
  {
    s=a+b;
    double v=s-a;
    e=(a-(s-v))+(b-v);
  }
  
  void three_summ(double &a,double &b,double &c,double x,double y,double z)
  {
    a=x;
    b=y;
    c=z;
    double t1,t2,t3;
    two_summ(t1,t2,a,b);
    two_summ(a,t3,c,t1);
    two_summ(b,c,t2,t3);
  }
  void three_summ(double &a,double &b,double &c)
  {three_summ(a,b,c,a,b,c);}
  
  void three_summ2(double &a,double &b,double x,double y,double c)
  {
    a=x;
    b=y;
    double t1,t2,t3;
    two_summ(t1,t2,a,b);
    two_summ(a,t3,c,t1);
    b=t2+t3;
  }
  void three_summ2(double &a,double &b,double c)
  {three_summ2(a,b,a,b,c);}
  
  double quick_three_accum(double &a,double &b,double c)
  {
    double s;
    bool za,zb;
    
    two_summ(s,b,b,c);
    two_summ(s,a,a,s);
    
    za=(a!=0);
    zb=(b!=0);
    
    if(za&&zb) return s;
    
    if(!zb)
      {
	b=a;
	a=s;
      }
    else a=s;
    
    return 0;
  }
  
  void split(double &a_hi,double &a_lo,double a)
  {
    double t=134217729*a;
    a_hi=t-(t-a);
    a_lo=a-a_hi;
    
    //master_printf("  splitting: %16.16lg==%16.16lg\n",((double)((int)a)),a);
    //master_printf("  splitting: %lg={%lg,%lg}\n",a,a_hi,a_lo);
  }
  
  void two_prod(double &p,double &e,double a,double b)
  {
    p=a*b;
    double a_hi,a_lo;
    double b_hi,b_lo;
    
    split(a_hi,a_lo,a);
    split(b_hi,b_lo,b);
    
    e=((a_hi*b_hi-p)+a_hi*b_lo+a_lo*b_hi)+a_lo*b_lo;
  }
  
  void renormalize(float_256 b,float_256_unr a)
  {
    //master_printf("Renormalizing: {%lg,%lg,%lg,%lg,%lg}\n",a[0],a[1],a[2],a[3],a[4]);
    double s0=0,s1=0,s2=0,s3=0;
    double c0,c1,c2,c3,c4;
    quick_two_summ(s0,c4,a[3],a[4]);
    quick_two_summ(s0,c3,a[2],s0);
    quick_two_summ(s0,c2,a[1],s0);
    quick_two_summ(c0,c1,a[0],s0);
    
    s0=c0;
    s1=c1;
    
    quick_two_summ(s0,s1,c0,c1);
    if(s1!=0)
      {
	quick_two_summ(s1,s2,s1,c2);
	if(s2!=0)
	  {
	    quick_two_summ(s2,s3,s2,c3);
	    if(s3!=0) s3+=c4;
	    else s2+=c4;
	  }
	else
	  {
	    quick_two_summ(s1,s2,s1,c3);
	    if(s2!=0) quick_two_summ(s2,s3,s2,c4);
	    else quick_two_summ(s1,s2,s1,c4);
	  }
      }
    else
      {
	quick_two_summ(s0,s1,s0,c2);
	if(s1!=0)
	  {
	    quick_two_summ(s1,s2,s1,c3);
	    if(s2!=0) quick_two_summ(s2,s3,s2,c4);
	    else quick_two_summ(s1,s2,s1,c4);
	  }
	else
	  {
	    quick_two_summ(s0,s1,s0,c3);
	    if(s1!=0) quick_two_summ(s1,s2,s1,c4);
	    else quick_two_summ(s0,s1,s0,c4);
	  }
      }
    
    b[0]=s0;
    b[1]=s1;
    b[2]=s2;
    b[3]=s3;
  }
  void renormalize(float_256 b,double a0,double a1,double a2,double a3,double a4)
  {
    float_256_unr a={a0,a1,a2,a3,a4};
    renormalize(b,a);
  }
  
  void renormalize(float_256 b,double c0,double c1,double c2,double c3)
  {
    double s0=0,s1=0,s2=0,s3=0;
    
    quick_two_summ(s0,c3,c2,c3);
    quick_two_summ(s0,c2,c1,s0);
    quick_two_summ(c0,c1,c0,s0);
    
    s0=c0;
    s1=c1;
    if(s1!=0)
      {
	quick_two_summ(s1,s2,s1,c2);
	if(s2!=0) quick_two_summ(s2,s3,s2,c3);
	else quick_two_summ(s1,s2,s1,c3);
      }
    else
      {
	quick_two_summ(s0,s1,s0,c2);
	if(s1!=0) quick_two_summ(s1,s2,s1,c3);
	else quick_two_summ(s0,s1,s0,c3);
      }
    
    b[0]=s0;
    b[1]=s1;
    b[2]=s2;
    b[3]=s3;
  }
  
  void float_256_summ_double(float_256 c,float_256 a,double b)
  {
    float_256_unr t;
    double e;
    
    two_summ(t[0],e,a[0],b);
    two_summ(t[1],e,a[1],e);
    two_summ(t[2],e,a[2],e);
    two_summ(t[3],t[4],a[3],e);
    
    renormalize(c,t);
  }
  
  void float_256_summ(float_256 c,float_256 a,float_256 b)
  {
    double s0,s1,s2,s3;
    double t0,t1,t2,t3;
    
    double v0,v1,v2,v3;
    double u0,u1,u2,u3;
    double w0,w1,w2,w3;
    
    s0=a[0]+b[0];
    s1=a[1]+b[1];
    s2=a[2]+b[2];
    s3=a[3]+b[3];
    
    v0=s0-a[0];
    v1=s1-a[1];
    v2=s2-a[2];
    v3=s3-a[3];
    
    u0=s0-v0;
    u1=s1-v1;
    u2=s2-v2;
    u3=s3-v3;
    
    w0=a[0]-u0;
    w1=a[1]-u1;
    w2=a[2]-u2;
    w3=a[3]-u3;
    
    u0=b[0]-v0;
    u1=b[1]-v1;
    u2=b[2]-v2;
    u3=b[3]-v3;
    
    t0=w0+u0;
    t1=w1+u1;
    t2=w2+u2;
    t3=w3+u3;
    
    two_summ(s1,t0,s1,t0);
    three_summ(s2,t1,t0);
    three_summ2(s3,t2,t1);
    t0=t0+t2+t3;
    
    renormalize(c,s0,s1,s2,s3,t0);
  }
  
  void float_256_summ_ieee(float_256 c,float_256 a,float_256 b)
  {
    int i=0,j=0,k=0;
    double s,t;
    double u,v;
    double x[4]={0,0,0,0};
    
    if(fabs(a[i])>fabs(b[j])) u=a[i++];
    else u=b[j++];
    if(fabs(a[i])>fabs(b[j])) v=a[i++];
    else v=b[j++];
    
    quick_two_summ(u,v,u,v);
    
    while(k<4)
      {
	if(i>=4&&j>=4)
	  {
	    x[k]=u;
	    if(k<3) x[++k]=v;
	    break;
	  }
	if(i>=4) t=b[j++];
	else
	  if(j>=4) t=a[i++];
	  else
	    if(fabs(a[i])>fabs(b[j])) t=a[i++];
	    else t=b[j++];
	
	s=quick_three_accum(u,v,t);
	
	if(s!=0) x[k++]=s;
      }
    
    for(k=i;k<4;k++) x[3]+=a[k];
    for(k=j;k<4;k++) x[3]+=b[k];
    
    renormalize(c,x[0],x[1],x[2],x[3]);
  }
  
  void float_256_prod(float_256 c,float_256 a,float_256 b)
  {
    double p0,p1,p2,p3,p4,p5;
    double q0,q1,q2,q3,q4,q5;
    double p6,p7,p8,p9;
    double q6,q7,q8,q9;
    double r0,r1;
    double t0,t1;
    double s0,s1,s2;
    
    two_prod(p0,q0,a[0],b[0]);
    two_prod(p1,q1,a[0],b[1]);
    two_prod(p2,q2,a[1],b[0]);
    two_prod(p3,q3,a[0],b[2]);
    two_prod(p4,q4,a[1],b[1]);
    two_prod(p5,q5,a[2],b[0]);
    
    three_summ(p1,p2,q0);
    
    three_summ(p2,q1,q2);
    three_summ(p3,p4,p5);
    two_summ(s0,t0,p2,p3);
    two_summ(s1,t1,q1,p4);
    s2=q2+p5;
    two_summ(s1,t0,s1,t0);
    s2+=(t0+t1);
    
    two_prod(p6,q6,a[0],b[3]);
    two_prod(p7,q7,a[1],b[2]);
    two_prod(p8,q8,a[2],b[1]);
    two_prod(p9,q9,a[3],b[0]);
    
    two_summ(q0,q3,q0,q3);
    two_summ(q4,q5,q4,q5);
    two_summ(p6,p7,p6,p7);
    two_summ(p8,p9,p8,p9);
    two_summ(t0,t1,q0,q4);
    t1+=(q3+q5);
    two_summ(r0,r1,p6,p8);
    r1+=(p7+p9);
    two_summ(q3,q4,t0,r0);
    q4+=(t1+r1);
    two_summ(t0,t1,q3,s1);
    t1+=q4;
    
    t1+=a[1]*b[3]+a[2]*b[2]+a[3]*b[1]+q6+q7+q8+q9+s2;
    
    renormalize(c,p0,p1,s0,t0,t1);
  }
  
  void float_256_prod_double(float_256 c,float_256 a,double b)
  {
    //master_printf("  taking product: {%lg %lg %lg %lg}*%lg\n",a[0],a[1],a[2],a[3],b);

    float_256_unr s;
    double p1,p2,p3;
    double q0,q1,q2;
    double t1,t2,t3;
    double r2;
    
    two_prod(s[0],q0,a[0],b);
    //master_printf("  a[0]: %16.16lg\n",a[0]);
    //master_printf("  b: %16.16lg\n",b);
    //master_printf("  s0: %lg\n",s[0]);
    //master_printf("  q0: %lg\n",q0);
    two_prod(p1,q1,a[1],b);
    //master_printf("  p1: %lg\n",p1);
    //master_printf("  q1: %lg\n",q1);
    two_prod(p2,q2,a[2],b);
    //master_printf("  p2: %lg\n",p2);
    //master_printf("  q2: %lg\n",q2);
    p3=a[3]*b;
    
    two_summ(s[1],t1,q0,p1);
    
    three_summ(s[2],t2,r2,p2,q1,t1);
    
    three_summ2(s[3],t3,p3,q2,t2);
    s[4]=r2+t3;
    
    //master_printf("  unrenormalized: %lg %lg %lg %lg %lg\n",s[0],s[1],s[2],s[3],s[4]);
    renormalize(c,s);
    //master_printf("  result: %lg %lg %lg %lg\n",c[0],c[1],c[2],c[3]);
  }
  
  void float_256_summassign(float_256 b,float_256 a)
  {float_256_summ(b,b,a);}
  void float_256_subt(float_256 c,float_256 a,float_256 b)
  {
    float_256 d;
    float_256_uminus(d,b);
    float_256_summ(c,a,d);
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
  {
    float_256 temp;
    float_256_prod(temp,out,in);
    float_256_copy(out,temp);
  }
  
  void float_256_div(float_256 c,float_256 a,float_256 b)
  {
    //master_printf("  a: %lg %lg %lg %lg\n",a[0],a[1],a[2],a[3]);
    //master_printf("  b: %lg %lg %lg %lg\n",b[0],b[1],b[2],b[3]);
    
    float_256 r;
    float_256_copy(r,a);
    
    float_256_unr q;
    
    for(int i=0;i<4;i++)
      {
	q[i]=r[0]/b[0];
	float_256 s;
	float_256_prod_double(s,b,q[i]);
	//master_printf("  subtracting: %lg %lg %lg %lg\n",s[0],s[1],s[2],s[3]);
	float_256_subtassign(r,s);
	//master_printf("  remaining: %lg %lg %lg %lg\n",r[0],r[1],r[2],r[3]);
      }
    
    q[4]=r[0]/b[0];
    //master_printf("  quotient: %lg %lg %lg %lg %lg\n",q[0],q[1],q[2],q[3],q[4]);
    renormalize(c,q);
    //float_256 s;
    //float_256_prod(s,b,c);
    //float_256_subt(r,a,s);
    //master_printf("  remnant: %lg %lg %lg %lg %lg\n",r[0],r[1],r[2],r[3]);
  }
  
  //integer power
  void float_256_pow_int(float_256 out,float_256 in,int d)
  {
    //master_printf("Taking: {%lg,%lg,%lg,%lg}^%d\n",in[0],in[1],in[2],in[3],d);
    
    //negative or null case
    if(d<=0)
      {
	float_256_from_double(out,1);
	//if negative
	if(d<0)
	  {
	    //compute inv and assign to out
	    float_256 inv;
	    //master_printf("Taking the inverse\n");
	    float_256_div(inv,out,in);
	    //master_printf("Copying\n");
	    float_256_copy(out,inv);
	    //multiply the remaining numer of times
	    for(int i=2;i<=-d;i++)
	      {
		//master_printf("Multiplying %d\n",i);
		float_256_prodassign(out,inv);
	      }
	  }
      }
    //positive case
    else
      {
	float_256_copy(out,in);
	for(int i=2;i<=d;i++) float_256_prodassign(out,in);
      }
  }
  
  //frac power
  void float_256_pow_int_frac(float_256 out,float_256 ext_in,int n,int d)
  {
    //master_printf("Computing (%16.16lg,%16.16lg,%16.16lg,%16.16lg)^(%d/%d)\n",ext_in[0],ext_in[1],ext_in[2],ext_in[3],n,d);
    
    //take a copy
    float_256 in;
    float_256_copy(in,ext_in);
    
    //compute by solving out^d=in^n=ref
    float_256 ref;
    float_256_pow_int(ref,in,n);
    //master_printf("Reference: (%lg,%lg,%lg,%lg)\n",ref[0],ref[1],ref[2],ref[3]);
    
    //let's start from a reasonable approx
    float_256_from_double(out,pow(in[0],(double)n/d));
    //another possible approach
    //float_256_from_128(out,1+(in[0]-1)*r);
    
    //comput the reciprocal of d
    float_256 recd,tempd;
    float_256_from_double(recd,1);
    float_256_from_double(tempd,d);
    float_256_div(recd,recd,tempd);
    //master_printf(" d=%d, 1/d: (%lg,%lg,%lg,%lg)\n",d,recd[0],recd[1],recd[2],recd[3]);

    //(out+err)^d=in^n -> err=out*rel_err, rel_err=(ref/out^d-1)/d
    int iter=0;
    float_256 rel_err;
    do
      {
	//compute out^d
	float_256 outd;
	float_256_pow_int(outd,out,d);
	//master_printf("Current guess: (%lg,%lg,%lg,%lg), pointing to: (%lg,%lg,%lg,%lg)\n",out[0],out[1],out[2],out[3],outd[0],outd[1],outd[2],outd[3]);
	
	//compute relative error
	float_256 temp;
	float_256_div(temp,ref,outd);
	float_256_prod(temp,temp,recd);
	float_256_subt(rel_err,temp,recd);
	
	//total err
	float_256 err;
	float_256_prod(err,rel_err,out);
	
	float_256_summassign(out,err);
	//master_printf("Iter %d err: %lg\n",iter,fabs(rel_err[0]));
	
	iter++;
      }
    while(fabs(rel_err[0])>2.7e-65);
  }
  
  void float_256_abs(float_256 a,float_256 b)
  {
    float_256 zero;
    float_256_put_to_zero(zero);
    if(float_256_is_greater(b,zero)) float_256_copy(a,b);
    else float_256_uminus(a,b);
  }
}
