#pragma once

#include <math.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace std;

class jack
{
public:
  int njack;
  double *data;
  void create(int);
  jack(const jack&);
  explicit jack(){data=NULL;njack=0;}
  explicit jack(int);
  explicit jack(int,int*);
  explicit jack(int,double*);
  ~jack();
  
  double operator[](int);
  jack operator=(double);
  jack operator=(const jack&);
  
  double med();
  double err();
  
  void reallocate_if_necessary(int);
  void put(double*);
  jack append_to_binfile(const char*,...);
  jack write_to_binfile(const char*,...);
};

//creation and assignment

void jack::put(double *in)
{
  memcpy(data,in,sizeof(double)*(njack+1));
}

void jack::create(int n)
{
  njack=n;
  data=new double[njack+1];
}

void jack::reallocate_if_necessary(int nj)
{
  if(njack!=nj)
    {
      if(data!=NULL) delete[] data;
      create(nj);
    }
}

jack::~jack(){if(data!=NULL) delete[]data;}

jack::jack(const jack &in) : njack(in.njack)
{
  data=new double[njack+1];
  put(in.data);
}

jack jack::operator=(const jack &in)
{
  reallocate_if_necessary(in.njack);
  put(in.data);
  
  return *this;
}

jack::jack(int n)
{
  create(n);
}

jack::jack(int n,double *in)
{
  create(n);
  put(in);
}

jack::jack(int n,int *in)
{
  double temp[n];
  for(int ijack=0;ijack<n+1;ijack++) temp[ijack]=in[ijack];
  (*this)=jack(n,temp);
}

// access and error

double jack::operator[](int ijack)
{
  return data[ijack];
}

jack jack::operator=(double a)
{
  for(int ijack=0;ijack<njack+1;ijack++) data[ijack]=a;
  return *this;
}

double jack::err()
{
  double sx=0,s2x=0;
  
  for(int ij=0;ij<njack;ij++)
    {
      sx+=data[ij];
      s2x+=data[ij]*data[ij];
    }
  sx/=njack;
  s2x/=njack;
  s2x-=sx*sx;
  
  return sqrt(s2x*(njack-1)); 
}

double jack::med()
{
  return data[njack];
}

jack jack::append_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"aw");
  int nw=fwrite(data,sizeof(double)*(njack+1),1,fout);
  if(nw!=1)
    {
      cerr<<"Error appending to file "<<buffer<<endl;
      exit(1);
    }
  fclose(fout);
  
  return *this;
}

jack jack::write_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"w");
  int nw=fwrite(data,sizeof(double)*(njack+1),1,fout);
  if(nw!=1)
    {
      cerr<<"Error appending to file "<<buffer<<endl;
      exit(1);
    }
  fclose(fout);
  
  return *this;
}

ostream& operator<<(ostream &out,const jack &obj)
{
  out<<jack(obj).med()<<" "<<jack(obj).err();
  return out;
}

//functions

double unary_minus(const double a){return -a;}
double unary_plus(const double a){return a;}
double sqr(const double a){return a*a;}

double double_summ(const double a,const double b){return a+b;}
double double_subt(const double a,const double b){return a-b;}
double double_prod(const double a,const double b){return a*b;}
double double_frac(const double a,const double b){return a/b;}

jack single_operator(const jack &a,double (*fun)(const double))
{
  int njack=a.njack;
  jack c(njack);
  transform(a.data,a.data+njack+1,c.data,ptr_fun(fun));

  return c;
}

jack pair_operator(const jack &a,const jack &b,double (*fun)(const double,const double))
{
  int njack=a.njack;
  if(b.njack!=njack)
    {
      cerr<<"Error, unmatched njack: "<<njack<<" "<<b.njack<<"!"<<endl;
      exit(1);
    }
  jack c(njack);
  transform(a.data,a.data+njack+1,b.data,c.data,ptr_fun(fun));

  return c;
}

jack pair_operator(const jack &a,const double b,double (*fun)(const double,const double),int du)
{
  int njack=a.njack;
  jack c(njack);
  
  for(int ijack=0;ijack<njack+1;ijack++)
    if(du==2) c.data[ijack]=fun(a.data[ijack],b);
    else      c.data[ijack]=fun(b,a.data[ijack]);

  return c;
}

jack pair_operator(const jack &a,const double b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
jack pair_operator(const double a,const jack &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

jack operator+(const jack &a,const jack &b){return pair_operator(a,b,double_summ);}
jack operator-(const jack &a,const jack &b){return pair_operator(a,b,double_subt);}
jack operator*(const jack &a,const jack &b){return pair_operator(a,b,double_prod);}
jack operator/(const jack &a,const jack &b){return pair_operator(a,b,double_frac);}

jack operator+=(jack &a,const jack &b){return a=a+b;}
jack operator-=(jack &a,const jack &b){return a=a-b;}
jack operator*=(jack &a,const jack &b){return a=a*b;}
jack operator/=(jack &a,const jack &b){return a=a/b;}

/////////////

jack operator+(const jack &a,const double b){return pair_operator(a,b,double_summ);}
jack operator-(const jack &a,const double b){return pair_operator(a,b,double_subt);}
jack operator*(const jack &a,const double b){return pair_operator(a,b,double_prod);}
jack operator/(const jack &a,const double b){return pair_operator(a,b,double_frac);}

jack operator+=(jack &a,const double b){return a=a+b;}
jack operator-=(jack &a,const double b){return a=a-b;}
jack operator*=(jack &a,const double b){return a=a*b;}
jack operator/=(jack &a,const double b){return a=a/b;}

jack operator+(const double a,const jack &b){return pair_operator(a,b,double_summ);}
jack operator-(const double a,const jack &b){return pair_operator(a,b,double_subt);}
jack operator*(const double a,const jack &b){return pair_operator(a,b,double_prod);}
jack operator/(const double a,const jack &b){return pair_operator(a,b,double_frac);}

jack operator+(const jack &a){return single_operator(a,unary_plus);}
jack operator-(const jack &a){return single_operator(a,unary_minus);}

jack sin(const jack &a){return single_operator(a,sin);}
jack cos(const jack &a){return single_operator(a,cos);}
jack tan(const jack &a){return single_operator(a,tan);}
jack asin(const jack &a){return single_operator(a,asin);}
jack acos(const jack &a){return single_operator(a,acos);}
jack atan(const jack &a){return single_operator(a,atan);}

jack sinh(const jack &a){return single_operator(a,sinh);}
jack cosh(const jack &a){return single_operator(a,cosh);}
jack tanh(const jack &a){return single_operator(a,tanh);}
jack asinh(const jack &a){return single_operator(a,asinh);}
jack acosh(const jack &a){return single_operator(a,acosh);}
jack atanh(const jack &a){return single_operator(a,atanh);}

jack sqr(const jack &a){return single_operator(a,sqr);}
jack sqrt(const jack &a){return single_operator(a,sqrt);}
jack pow(const jack &a,double b){return pair_operator(a,b,pow);}
