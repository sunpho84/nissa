#pragma once

#include <math.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace std;

class jvec
{
public:
  int nel;
  int njack;
  jack *data;
  void create(int,int);
  jvec();
  jvec(const jvec&);
  explicit jvec(int,int);
  explicit jvec(int,double*);
  
  void put(double *);
  jvec load(const char *,int);
  jvec load_naz(const char *,int);
  
  jack& operator[](int);
  jvec operator=(double);
};

//creation and assignment
void jvec::create(int ne,int nj)
{
  nel=ne;
  njack=nj;
  data=new jack[nel];
  for(int iel=0;iel<nel;iel++) data[iel].create(nj);
}

jvec::jvec(){}

jvec::jvec(const jvec &in) : nel(in.nel),njack(in.njack)
{
  create(nel,njack);
  for(int iel=0;iel<nel;iel++) data[iel]=in.data[iel];
}

jvec::jvec(int ne,int nj)
{
  create(ne,nj);
}

void jvec::put(double *in)
{
  for(int iel=0;iel<nel;iel++)
    data[iel].put(in+iel*(njack+1));
}

jvec jvec::load(const char *path,int i)
{
  double in[nel*(njack+1)];
  
  FILE *fin=open_file(path,"r");

  if(fseek(fin,i*sizeof(double)*nel*(njack+1),SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",i);
      exit(1);
    }
  int stat=fread(in,sizeof(double)*nel*(njack+1),1,fin);
  if(stat!=1)
    {
      if(stat==EOF)
	{
	  fprintf(stderr,"Error, reached EOF while reading data!\n");
	  exit(1);
	}
      else
        {
	  perror("Error while reading data!");
	  exit(1);
	}
    }
  
  put(in);
  
  fclose(fin);
  
  return (*this);
}

jvec jvec::load_naz(const char *path,int icorr)
{
  double in[nel][2][njack+1];
  
  FILE *fin=open_file(path,"r");

  if(fseek(fin,2*(icorr/2)*sizeof(double)*nel*(njack+1),SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",icorr);
      exit(1);
    }
  int stat=fread(in,sizeof(double)*2*nel*(njack+1),1,fin);
  if(stat!=1)
    {
      if(stat==EOF)
	{
	  fprintf(stderr,"Error, reached EOF while reading data!\n");
	  exit(1);
	}
      else
        {
	  perror("Error while reading data!\n");
	  exit(1);
	}
    }
  
  int ri=icorr%2;
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      data[iel].data[ijack]=in[iel][ri][ijack];
  
  fclose(fin);
  
  return (*this);
}

jvec jvec_load(const char *path,int nel,int njack,int i)
{
  jvec out(nel,njack);
  
  return out.load(path,i);
}

jvec jvec_load_naz(const char *path,int nel,int njack,int i)
{
  jvec out(nel,njack);
  
  return out.load_naz(path,i);
}

ostream& operator<<(ostream &out,const jvec &obj)
{
  for(int iel=0;iel<obj.nel;iel++) out<<iel<<" "<<obj.data[iel]<<endl;
  
  return out;
}

jack& jvec::operator[](int iel)
{
  return data[iel];
}

jvec jvec::operator=(double in)
{
  for(int iel=0;iel<nel;iel++) data[iel]=0;
  return *this;
}

jvec single_operator(const jvec &a,double (*fun)(const double))
{
  int nel=a.nel;
  int njack=a.njack;
  jvec c(nel,njack);
  
  for(int iel=0;iel<nel;iel++) c.data[iel]=single_operator(a.data[iel],fun);

  return c;
}

jvec pair_operator(const jvec &a,const jvec &b,double (*fun)(const double,const double))
{
  int nel=a.nel;
  int njack=a.njack;
  if(b.njack!=njack||b.nel!=nel)
    {
      cerr<<"Error, unmatched njack or nel!"<<endl;
      exit(1);
    }
  jvec c(nel,njack);
  for(int iel=0;iel<nel;iel++)
    c.data[iel]=pair_operator(a.data[iel],b.data[iel],fun);
  
  return c;
}

jvec pair_operator(const jvec &a,const double b,double (*fun)(const double,const double),int du)
{
  int nel=a.nel;
  int njack=a.njack;
  jvec c(nel,njack);
  
  for(int iel=0;iel<nel;iel++) c.data[iel]=pair_operator(a.data[iel],b,fun,du);
  
  return c;
}

jvec pair_operator(const jvec &a,const double b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
jvec pair_operator(const double a,const jvec &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

jvec pair_operator(const jvec &a,const jack b,double (*fun)(const double,const double),int du)
{
  int nel=a.nel;
  int njack=a.njack;
  jvec c(nel,njack);
  
  for(int iel=0;iel<nel;iel++)
    if(du==2) c.data[iel]=pair_operator(a.data[iel],b,fun);
    else      c.data[iel]=pair_operator(b,a.data[iel],fun);
  
  return c;
}

jvec pair_operator(const jvec &a,const jack b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
jvec pair_operator(const jack a,const jvec &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

jvec operator+(const jvec &a,const jvec &b){return pair_operator(a,b,double_summ);}
jvec operator-(const jvec &a,const jvec &b){return pair_operator(a,b,double_subt);}
jvec operator*(const jvec &a,const jvec &b){return pair_operator(a,b,double_prod);}
jvec operator/(const jvec &a,const jvec &b){return pair_operator(a,b,double_frac);}

jvec operator+=(const jvec &a,const jvec &b){return a+b;}
jvec operator-=(const jvec &a,const jvec &b){return a-b;}
jvec operator*=(const jvec &a,const jvec &b){return a*b;}
jvec operator/=(const jvec &a,const jvec &b){return a/b;}

///////////

jvec operator+(const jvec &a,const jack &b){return pair_operator(a,b,double_summ);}
jvec operator-(const jvec &a,const jack &b){return pair_operator(a,b,double_subt);}
jvec operator*(const jvec &a,const jack &b){return pair_operator(a,b,double_prod);}
jvec operator/(const jvec &a,const jack &b){return pair_operator(a,b,double_frac);}

jvec operator+(const jack &a,const jvec &b){return pair_operator(a,b,double_summ);}
jvec operator-(const jack &a,const jvec &b){return pair_operator(a,b,double_subt);}
jvec operator*(const jack &a,const jvec &b){return pair_operator(a,b,double_prod);}
jvec operator/(const jack &a,const jvec &b){return pair_operator(a,b,double_frac);}

jvec operator+=(jvec &a,const jack &b){return a=a+b;}
jvec operator-=(jvec &a,const jack &b){return a=a-b;}
jvec operator*=(jvec &a,const jack &b){return a=a*b;}
jvec operator/=(jvec &a,const jack &b){return a=a/b;}

////////////

jvec operator+(const jvec &a,const double b){return pair_operator(a,b,double_summ);}
jvec operator-(const jvec &a,const double b){return pair_operator(a,b,double_subt);}
jvec operator*(const jvec &a,const double b){return pair_operator(a,b,double_prod);}
jvec operator/(const jvec &a,const double b){return pair_operator(a,b,double_frac);}

jvec operator+=(jvec &a,const double b){return a=a+b;}
jvec operator-=(jvec &a,const double b){return a=a-b;}
jvec operator*=(jvec &a,const double b){return a=a*b;}
jvec operator/=(jvec &a,const double b){return a=a/b;}

jvec operator+(const double a,const jvec &b){return pair_operator(a,b,double_summ);}
jvec operator-(const double a,const jvec &b){return pair_operator(a,b,double_subt);}
jvec operator*(const double a,const jvec &b){return pair_operator(a,b,double_prod);}
jvec operator/(const double a,const jvec &b){return pair_operator(a,b,double_frac);}

////////////

jvec operator+(const jvec &a){return single_operator(a,unary_plus);}
jvec operator-(const jvec &a){return single_operator(a,unary_minus);}
jvec sin(const jvec &a){return single_operator(a,sin);}
jvec cos(const jvec &a){return single_operator(a,cos);}
jvec tan(const jvec &a){return single_operator(a,tan);}
jvec asin(const jvec &a){return single_operator(a,asin);}
jvec acos(const jvec &a){return single_operator(a,acos);}
jvec atan(const jvec &a){return single_operator(a,atan);}

jvec sqrt(const jvec &a){return single_operator(a,sqrt);}
jvec pow(const jvec &a,double b){return pair_operator(a,b,pow);}

//////////////
