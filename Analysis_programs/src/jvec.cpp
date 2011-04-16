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
  explicit jvec(){data=NULL;nel=0;njack=0;}
  jvec(const jvec&);
  jvec(int,int);

  void reallocate_if_necessary(int,int);
  void put(double *);
  jvec load(const char *,int);
  jvec load_naz(const char *,int);
  void print_to_file(const char *,...);
  
  jack& operator[](int);
  jvec operator=(double in){for(int iel=0;iel<nel;iel++) data[iel]=in;return *this;}
  jvec operator=(const jvec&);

  jvec first_half();
  jvec simmetric();
  jvec simmetrized(int parity);
};

ostream& operator<<(ostream &out,const jvec &obj);

//creation and assignment
void jvec::create(int ne,int nj)
{
  nel=ne;
  njack=nj;
  data=new jack[nel];
  for(int iel=0;iel<nel;iel++) data[iel].create(nj);
}

void jvec::reallocate_if_necessary(int ne,int nj)
{
  if(nel!=ne||njack!=nj)
    {
      if(data!=NULL) delete[] data;
      create(ne,nj);
    }
}

jvec::jvec(const jvec &in) : nel(in.nel),njack(in.njack)
{
  create(nel,njack);
  for(int iel=0;iel<nel;iel++) data[iel]=in.data[iel];
}

jvec jvec::operator=(const jvec &in)
{
  reallocate_if_necessary(in.nel,in.njack);
  for(int iel=0;iel<nel;iel++) data[iel]=in.data[iel];
  return *this;
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

void jvec::print_to_file(const char *format,...)
{
  char buffer[1024];
  va_list args;

  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);

  ofstream fout(buffer);
  fout<<(*this);
  fout.close();
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

jvec operator+=(jvec &a,const jvec &b){return a=a+b;}
jvec operator-=(jvec &a,const jvec &b){return a=a-b;}
jvec operator*=(jvec &a,const jvec &b){return a=a*b;}
jvec operator/=(jvec &a,const jvec &b){return a=a/b;}

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

jvec sinh(const jvec &a){return single_operator(a,sinh);}
jvec cosh(const jvec &a){return single_operator(a,cosh);}
jvec tanh(const jvec &a){return single_operator(a,tanh);}
jvec asinh(const jvec &a){return single_operator(a,asinh);}
jvec acosh(const jvec &a){return single_operator(a,acosh);}
jvec atanh(const jvec &a){return single_operator(a,atanh);}

jvec sqrt(const jvec &a){return single_operator(a,sqrt);}
jvec pow(const jvec &a,double b){return pair_operator(a,b,pow);}

//////////////

jvec jvec::first_half()
{
  if(nel%2!=0)
    {
      cerr<<"Error, required the first half of an odd-length ("<<nel<<") jvector!"<<endl;
      exit(1);
    }
  
  jvec c(nel/2,njack);
  for(int iel=0;iel<nel/2;iel++)
    c.data[iel]=data[iel];
  
  return c;
}

jvec jvec::simmetric()
{
  jvec c(nel,njack);
  for(int iel=0;iel<nel;iel++)
    c.data[iel]=data[nel-iel-1];

  return c;
}

jvec jvec::simmetrized(int parity)
{
  if(abs(parity)!=1)
    {
      cerr<<"Error, parity required for simmetrization: "<<parity<<endl;
      exit(1);
    }
  
  if(nel%2!=0)
    {
      cerr<<"Error, required to simmetrize an odd-length jvector!"<<endl;
      exit(1);
    }
  
  jvec c(nel/2+1,njack);
  
  for(int iel=0;iel<=nel/2;iel++)
    {
      if(parity==1) c.data[iel]=(data[iel]+data[(nel-iel)%nel])/2;
      else          c.data[iel]=(data[iel]-data[(nel-iel)%nel])/2;
    }

  return c;
}

jack constant_fit(jvec in,int tin,int tfin)
{
  int njack=in.njack;
  jack E(njack),temp;

  E=0;
  double norm=0;
  for(int iel=tin;iel<=tfin;iel++)
    {
      double err=in.data[iel].err();
      double weight=1/(err*err);
      E+=in[iel]*weight;
      temp=in[iel]*weight;
      norm+=weight;
    }
  return E/norm;
}

void linear_fit(jvec in,jack &m,jack &q,int tin,int tfin)
{
  int njack=in.njack;
  double S,Sx,Sx2;
  jack Sxy(njack),Sy(njack);
  
  Sx2=S=Sx=0;
  Sxy=Sy=0;
  for(int iel=tin;iel<=tfin;iel++)
    {
      double err=in.data[iel].err();
      double weight=1/(err*err);
      int x=iel;
      jack y=in.data[iel];

      S+=weight;
      Sx+=x*weight;
      Sx2+=x*x*weight;
      Sxy+=x*y*weight;
      Sy+=y*weight;
    }
  
  double delta=S*Sx2-Sx*Sx;
  m=(S*Sxy-Sx*Sy)/delta;
  q=(Sx2*Sy-Sxy*Sx)/delta;
}
