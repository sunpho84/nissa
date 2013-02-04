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

void three_summ(double &a,double &b,double &c)
{
  double t1,t2,t3;
  two_summ(t1,t2,a,b);
  two_summ(a,t3,c,t1);
  two_summ(b,c,t2,t3);
}

void three_summ2(double &a,double &b,double &c)
{
  double t1,t2,t3;
  two_summ(t1,t2,a,b);
  two_summ(a,t3,c,t1);
  b=t2+t3;
}

void split(double &a_hi,double &a_lo,double a)
{
  double t=134217729*a;
  a_hi=t-(t-a);
  a_lo=a-a_hi;
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
  //printf("before: ");
  //for(int i=0;i<5;i++) printf("%lg ",a[i]);
  //printf("\n");
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
  
  //printf("after: ");
  //for(int i=0;i<4;i++) printf("%lg ",b[i]);
  //printf("\n");

}
void renormalize(float_256 b,double a0,double a1,double a2,double a3,double a4)
{
  float_256_unr a={a0,a1,a2,a3,a4};
  renormalize(b,a);
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
  three_summ(s2,t0,t1);
  three_summ2(s3,t0,t2);
  t0=t0+t1+t3;

  renormalize(c,s0,s1,s2,s3,t0);
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
  float_256_unr s;
  double p0,p1,p2,p3;
  double q0,q1,q2;

  two_prod(p0,q0,a[0],b);
  two_prod(p1,q1,a[1],b);
  two_prod(p2,q2,a[2],b);
  p3=a[3]*b;
  
  s[0]=p0;
  
  two_summ(s[1],s[2],q0,p1);
  
  three_summ(s[2],q1,p2);
  
  three_summ2(q1,q2,p3);
  s[3]=q1;  
  s[4]=q2+p2;
  
  renormalize(c,s);
}
