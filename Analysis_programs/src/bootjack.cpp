#include <math.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace std;

class TYPE
{
public:
  int njack;
#ifdef BOOT
  int nboot;
#endif
  double *data;

#ifdef BOOT
  void create(int a,int b){nboot=a;njack=b;data=new double[a+1];}
  TYPE(const TYPE& in){create(in.nboot,in.njack);put(in.data);}
  explicit TYPE(){data=NULL;nboot=njack=0;}
  explicit TYPE(int a,int b){create(a,b);}
  explicit TYPE(int nb,int nj,int *in){double temp[nb+1];for(int i=0;i<nb+1;i++) temp[i]=in[i];TYPE(nb,nj,temp);}
  explicit TYPE(int nb,int nj,double *in){create(nb,nj);put(in);}
  void reallocate_if_necessary(int nb,int nj){if(nboot!=nb||njack!=nj){if(data!=NULL) delete[] data;create(nb,nj);}}
  TYPE operator=(const TYPE& in){reallocate_if_necessary(in.nboot,in.njack);put(in.data);return *this;}
#else
  void create(int nj){njack=nj;data=new double[nj+1];}
  TYPE(const TYPE& in){create(in.njack);put(in.data);}
  explicit TYPE(){data=NULL;njack=0;}
  explicit TYPE(int nj){create(nj);}
  explicit TYPE(int nj,int *in){double temp[nj+1];for(int i=0;i<nj+1;i++) temp[i]=in[i];TYPE(nj,temp);}
  explicit TYPE(int nj,double *in){create(nj);put(in);}
  void reallocate_if_necessary(int nj){if(njack!=nj){if(data!=NULL) delete[] data;create(nj);}}
  TYPE operator=(const TYPE& in){reallocate_if_necessary(in.njack);put(in.data);return *this;}
#endif
  ~TYPE(){if(data!=NULL) delete[]data;data=NULL;}
  
  double operator[](int i){return data[i];}
  TYPE operator=(double a){for(int i=0;i<N+1;i++) data[i]=a;return *this;}
  
  double med(){return data[N];}
  double err();
  
  void fill_gauss(double med,double sigma){for(int itype=0;itype<N;itype++)data[itype]=rng(med,sigma/sqrt(njack-1));data[N]=med;}
  void put(double* in){memcpy(data,in,sizeof(double)*(N+1));}

  TYPE append_to_binfile(const char*,...);
  TYPE write_to_binfile(const char*,...);
};

double TYPE::err()
{
  double sx=0,s2x=0;
  
  for(int ij=0;ij<N;ij++)
    {
      sx+=data[ij];
      s2x+=data[ij]*data[ij];
    }
  sx/=N;
  s2x/=N;
  s2x-=sx*sx;
  
  return sqrt(s2x*(njack-1));
}

TYPE TYPE::append_to_binfile(const char *format,...)
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

TYPE TYPE::write_to_binfile(const char *format,...)
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

ostream& operator<<(ostream &out,const TYPE &obj)
{
  out<< TYPE(obj).med()<<" "<< TYPE(obj).err();
  return out;
}

//functions

TYPE single_operator(const TYPE &a,double (*fun)(const double))
{
  int N=a.N;
  TYPE c(a);
  transform(a.data,a.data+N+1,c.data,ptr_fun(fun));

  return c;
}

TYPE pair_operator(const TYPE &a,const TYPE &b,double (*fun)(const double,const double))
{
  int N=a.N;
  if(b.N!=N)
    {
      cerr<<"Error, unmatched N: "<<N<<" "<<b.N<<"!"<<endl;
      exit(1);
    }
  TYPE c(a);
  transform(a.data,a.data+N+1,b.data,c.data,ptr_fun(fun));

  return c;
}

TYPE pair_operator(const TYPE &a,const double b,double (*fun)(const double,const double),int du)
{
  int N=a.N;
  TYPE c(a);
  
  for(int i=0;i<N+1;i++)
    if(du==2) c.data[i]=fun(a.data[i],b);
    else      c.data[i]=fun(b,a.data[i]);
  
  return c;
}

TYPE pair_operator(const TYPE &a,const double b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
TYPE pair_operator(const double a,const TYPE &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

TYPE operator+(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_summ);}
TYPE operator-(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_subt);}
TYPE operator*(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_prod);}
TYPE operator/(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_frac);}

TYPE operator+=(TYPE &a,const TYPE &b){return a=a+b;}
TYPE operator-=(TYPE &a,const TYPE &b){return a=a-b;}
TYPE operator*=(TYPE &a,const TYPE &b){return a=a*b;}
TYPE operator/=(TYPE &a,const TYPE &b){return a=a/b;}

/////////////

TYPE operator+(const TYPE &a,const double b){return pair_operator(a,b,double_summ);}
TYPE operator-(const TYPE &a,const double b){return pair_operator(a,b,double_subt);}
TYPE operator*(const TYPE &a,const double b){return pair_operator(a,b,double_prod);}
TYPE operator/(const TYPE &a,const double b){return pair_operator(a,b,double_frac);}

TYPE operator+=( TYPE &a,const double b){return a=a+b;}
TYPE operator-=( TYPE &a,const double b){return a=a-b;}
TYPE operator*=( TYPE &a,const double b){return a=a*b;}
TYPE operator/=( TYPE &a,const double b){return a=a/b;}

TYPE operator+(const double a,const TYPE &b){return pair_operator(a,b,double_summ);}
TYPE operator-(const double a,const TYPE &b){return pair_operator(a,b,double_subt);}
TYPE operator*(const double a,const TYPE &b){return pair_operator(a,b,double_prod);}
TYPE operator/(const double a,const TYPE &b){return pair_operator(a,b,double_frac);}

TYPE operator+(const TYPE &a){return single_operator(a,unary_plus);}
TYPE operator-(const TYPE &a){return single_operator(a,unary_minus);}

TYPE sin(const TYPE &a){return single_operator(a,sin);}
TYPE cos(const TYPE &a){return single_operator(a,cos);}
TYPE tan(const TYPE &a){return single_operator(a,tan);}
TYPE asin(const TYPE &a){return single_operator(a,asin);}
TYPE acos(const TYPE &a){return single_operator(a,acos);}
TYPE atan(const TYPE &a){return single_operator(a,atan);}

TYPE sinh(const TYPE &a){return single_operator(a,sinh);}
TYPE cosh(const TYPE &a){return single_operator(a,cosh);}
TYPE tanh(const TYPE &a){return single_operator(a,tanh);}
TYPE asinh(const TYPE &a){return single_operator(a,asinh);}
TYPE acosh(const TYPE &a){return single_operator(a,acosh);}
TYPE atanh(const TYPE &a){return single_operator(a,atanh);}

TYPE sqr(const TYPE &a){return single_operator(a,sqr);}
TYPE sqrt(const TYPE &a){return single_operator(a,sqrt);}
TYPE pow(const TYPE &a,double b){return pair_operator(a,b,pow);}
