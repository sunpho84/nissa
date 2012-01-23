#include <math.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>

using namespace std;

class VTYPE
{
public:
  int nel;
#ifdef BVEC
  int nboot;
#endif
  int njack;
  TYPE *data;
  
#ifdef BVEC
  void create(int ne,int nb,int nj){nel=ne;nboot=nb;njack=nj;data=new boot[nel];for(int iel=0;iel<nel;iel++) data[iel].create(nb,nj);}
  VTYPE(){data=NULL;nel=nboot=njack=0;}
  VTYPE(const VTYPE& a) : nel(a.nel),nboot(a.nboot),njack(a.njack){create(nel,nboot,njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];}
  VTYPE(int ne,int nb,int nj){create(ne,nb,nj);}
  void reallocate_if_necessary(int ne,int nb,int nj){if(nel!=ne||njack!=nj||nb!=nboot){delete[] data;create(ne,nb,nj);}}
  VTYPE operator=(const VTYPE& a){
    reallocate_if_necessary(a.nel,a.nboot,a.njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];return *this;}
#else
  void create(int ne,int nj){nel=ne;njack=nj;data=new jack[nel];for(int iel=0;iel<nel;iel++) data[iel].create(nj);}
  explicit VTYPE(){data=NULL;nel=njack=0;}
  VTYPE(const VTYPE& a) : nel(a.nel),njack(a.njack){create(nel,njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];}
  VTYPE(int ne,int nj){create(ne,nj);}
  void reallocate_if_necessary(int ne,int nj){if(nel!=ne||njack!=nj){if(data!=NULL)delete[] data;create(ne,nj);}}
  VTYPE operator=(const VTYPE& a){
    reallocate_if_necessary(a.nel,a.njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];return *this;}
#endif
  ~VTYPE(){if(data!=NULL) delete[] data;data=NULL;}
  
  void put(double **out){for(int iel=0;iel<nel;iel++) data[iel].put(out[iel]);}
  void put(double *out){for(int iel=0;iel<nel;iel++) data[iel].put(out+iel*(N+1));}
  void get(double *in){for(int iel=0;iel<nel;iel++) data[iel].get(in+iel*(N+1));}
  VTYPE load(const char *,int);
  VTYPE load(FILE *,int);
  VTYPE load_naz(const char *,int);
  void print_to_file(const char *,...);
  
  TYPE& operator[](int i){return data[i];}
  VTYPE operator=(double in){for(int iel=0;iel<nel;iel++) data[iel]=in;return *this;}
  VTYPE operator=(const TYPE& a){
    if(data==NULL){cerr<<"Error, using unallocated vector!"<<endl;exit(1);}for(int iel=0;iel<nel;iel++) data[iel]=a;return *this;}

  
  VTYPE first_half();
  VTYPE simmetric();
  VTYPE simmetrized(int parity);
  
  VTYPE write_to_binfile(FILE *);
  VTYPE append_to_binfile(const char*,...);
  VTYPE write_to_binfile(const char*,...);
};

ostream& operator<<(ostream &out,const VTYPE &obj);

VTYPE VTYPE::write_to_binfile(FILE *fout)
{
  double out[nel*(N+1)];
  get(out);
  int nw=fwrite(out,sizeof(double),(N+1)*nel,fout);
  if(nw!=(N+1)*nel)
    {
      cerr<<"Error writing to file"<<endl;
      exit(1);
    }
  
  return *this;
}

VTYPE VTYPE::append_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"aw");
  write_to_binfile(fout);
  fclose(fout);
  
  return *this;
}

VTYPE VTYPE::write_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"w");
  write_to_binfile(fout);
  fclose(fout);
  
  return *this;
}

VTYPE VTYPE::load(FILE *fin,int i)
{
  double in[nel*(N+1)];
  
  if(fseeko(fin,(off_t)i*sizeof(double)*nel*(N+1),SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",i);
      exit(1);
    }
  
  int stat=fread(in,sizeof(double),nel*(N+1),fin);
  if(stat!=nel*(N+1))
    {
      if(stat==EOF)
	{
	  fprintf(stderr,"Error, reached EOF while reading data!\n");
	  exit(1);
	}
      else
        {
	  fprintf(stderr,"Error while reading data from file, obtained: %d while reading %d elements\n",stat,nel*(N+1));
	  exit(1);
	}
    }
  
  put(in);
  
  return *this;
}  

VTYPE VTYPE::load(const char *path,int i)
{
  cout<<"Loading from path "<<path<<endl;
  FILE *fin=open_file(path,"r");
  
  load(fin,i);
  fclose(fin);
  
  return (*this);
}

void VTYPE::print_to_file(const char *format,...)
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

VTYPE VTYPE::load_naz(const char *path,int icorr)
{
  double in[nel][2][N+1];
  
  FILE *fin=open_file(path,"r");
  
  int off=(off_t)2*(icorr/2)*sizeof(double)*nel*(N+1);
  if(fseeko(fin,off,SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",icorr);
      exit(1);
    }
  int stat=fread(in,sizeof(double)*2*nel*(N+1),1,fin);
  if(stat!=1)
    {
      if(stat==EOF)
	{
	  fprintf(stderr,"Error, reached EOF while reading data!\n");
	  exit(1);
	}
      else
        {
	  perror("Error while reading data");
	  exit(1);
	}
    }
  
  int ri=icorr%2;
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<N+1;ijack++)
      data[iel].data[ijack]=in[iel][ri][ijack];
  
  fclose(fin);
  
  return (*this);
}

#ifdef BVEC
bvec bvec_load(const char *path,int nel,int nboot,int njack,int i)
{
  bvec out(nel,nboot,njack);
  
  return out.load(path,i);
}

bvec bvec_load_naz(const char *path,int nel,int nboot,int njack,int i)
{
  bvec out(nel,nboot,njack);
  
  return out.load_naz(path,i);
}
#else
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
#endif

ostream& operator<<(ostream &out,const VTYPE &obj)
{
  for(int iel=0;iel<obj.nel;iel++) 
    {
      double med=obj.data[iel].med();
      double err=obj.data[iel].err();
      if(!isnan(med) && !isnan(err))
	out<<iel<<" "<<obj.data[iel]<<endl;
    }
  
  return out;
}

VTYPE single_operator(const VTYPE &a,double (*fun)(const double))
{
  int nel=a.nel;
  int njack=a.njack;
#ifdef BVEC
  int nboot=a.nboot;
  VTYPE c(nel,nboot,njack);
#else
  VTYPE c(nel,njack);
#endif

  for(int iel=0;iel<nel;iel++) c.data[iel]=single_operator(a.data[iel],fun);

  return c;
}

VTYPE pair_operator(const VTYPE &a,const VTYPE &b,double (*fun)(const double,const double))
{
  int nel=a.nel;
  int njack=a.njack;
#ifdef BVEC
  int nboot=a.nboot;
  if(b.nboot!=nboot){cerr<<"Error, unmatched nboot!"<<endl;exit(1);}
  if(b.njack!=njack){cerr<<"Error, unmatched njack!"<<endl;exit(1);}
  if(b.nel!=nel){cerr<<"Error, unmatched nel!"<<endl;exit(1);}
  bvec c(nel,nboot,njack);
#else
  if(b.njack!=njack||b.nel!=nel){cerr<<"Error, unmatched njack or nel!"<<endl;exit(1);}
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++) c.data[iel]=pair_operator(a.data[iel],b.data[iel],fun);
  
  return c;
}

VTYPE pair_operator(const VTYPE &a,const double b,double (*fun)(const double,const double),int du)
{
  int nel=a.nel;
  int njack=a.njack;
#ifdef BVEC
  int nboot=a.nboot;
  bvec c(nel,nboot,njack);
#else
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++) c.data[iel]=pair_operator(a.data[iel],b,fun,du);
  
  return c;
}

VTYPE pair_operator(const VTYPE &a,const double b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
VTYPE pair_operator(const double a,const VTYPE &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

VTYPE pair_operator(const VTYPE &a,const TYPE b,double (*fun)(const double,const double),int du)
{
  int nel=a.nel;
#ifdef BVEC
  bvec c(nel,a.nboot,a.njack);
#else
  jvec c(nel,a.njack);
#endif
  
  for(int iel=0;iel<nel;iel++)
    if(du==2) c.data[iel]=pair_operator(a.data[iel],b,fun);
    else      c.data[iel]=pair_operator(b,a.data[iel],fun);
  
  return c;
}

VTYPE pair_operator(const VTYPE &a,const TYPE b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
VTYPE pair_operator(const TYPE a,const VTYPE &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

VTYPE operator+(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_frac);}

VTYPE operator+=(VTYPE &a,const VTYPE &b){return a=a+b;}
VTYPE operator-=(VTYPE &a,const VTYPE &b){return a=a-b;}
VTYPE operator*=(VTYPE &a,const VTYPE &b){return a=a*b;}
VTYPE operator/=(VTYPE &a,const VTYPE &b){return a=a/b;}

///////////

VTYPE operator+(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_frac);}

VTYPE operator+(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_frac);}

VTYPE operator+=(VTYPE &a,const TYPE &b){return a=a+b;}
VTYPE operator-=(VTYPE &a,const TYPE &b){return a=a-b;}
VTYPE operator*=(VTYPE &a,const TYPE &b){return a=a*b;}
VTYPE operator/=(VTYPE &a,const TYPE &b){return a=a/b;}

////////////

VTYPE operator+(const VTYPE &a,const double b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const VTYPE &a,const double b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const VTYPE &a,const double b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const VTYPE &a,const double b){return pair_operator(a,b,double_frac);}

VTYPE operator+=(VTYPE &a,const double b){return a=a+b;}
VTYPE operator-=(VTYPE &a,const double b){return a=a-b;}
VTYPE operator*=(VTYPE &a,const double b){return a=a*b;}
VTYPE operator/=(VTYPE &a,const double b){return a=a/b;}

VTYPE operator+(const double a,const VTYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const double a,const VTYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const double a,const VTYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const double a,const VTYPE &b){return pair_operator(a,b,double_frac);}

////////////

VTYPE operator+(const VTYPE &a){return single_operator(a,unary_plus);}
VTYPE operator-(const VTYPE &a){return single_operator(a,unary_minus);}

VTYPE sin(const VTYPE &a){return single_operator(a,sin);}
VTYPE cos(const VTYPE &a){return single_operator(a,cos);}
VTYPE tan(const VTYPE &a){return single_operator(a,tan);}
VTYPE asin(const VTYPE &a){return single_operator(a,asin);}
VTYPE acos(const VTYPE &a){return single_operator(a,acos);}
VTYPE atan(const VTYPE &a){return single_operator(a,atan);}

VTYPE sinh(const VTYPE &a){return single_operator(a,sinh);}
VTYPE cosh(const VTYPE &a){return single_operator(a,cosh);}
VTYPE tanh(const VTYPE &a){return single_operator(a,tanh);}
VTYPE asinh(const VTYPE &a){return single_operator(a,asinh);}
VTYPE acosh(const VTYPE &a){return single_operator(a,acosh);}
VTYPE atanh(const VTYPE &a){return single_operator(a,atanh);}

VTYPE sqr(const VTYPE &a){return single_operator(a,sqr);}
VTYPE log(const VTYPE &a){return single_operator(a,log);}
VTYPE sqrt(const VTYPE &a){return single_operator(a,sqrt);}
VTYPE pow(const VTYPE &a,double b){return pair_operator(a,b,pow);}

//////////////

VTYPE VTYPE::first_half()
{
  if(nel%2!=0)
    {
      cerr<<"Error, required the first half of an odd-length ("<<nel<<") VTYPEtor!"<<endl;
      exit(1);
    }
  
#ifdef BVEC
  bvec c(nel/2+1,nboot,njack);
#else
  jvec c(nel/2+1,njack);
#endif
  for(int iel=0;iel<=nel/2;iel++)
    c.data[iel]=data[iel];
  
  return c;
}

VTYPE VTYPE::simmetric()
{
#ifdef BVEC
  bvec c(nel,nboot,njack);
#else
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++)
    c.data[iel]=data[(nel-iel)%nel];

  return c;
}

VTYPE VTYPE::simmetrized(int parity)
{
  if(abs(parity)!=1)
    {
      cerr<<"Error, parity required for simmetrization: "<<parity<<endl;
      exit(1);
    }
  
  if(nel%2!=0)
    {
      cerr<<"Error, required to simmetrize an odd-length VTYPEtor!"<<endl;
      exit(1);
    }
  
#ifdef BVEC
  bvec c(nel/2+1,nboot,njack);
#else
  jvec c(nel/2+1,njack);
#endif
  
  for(int iel=0;iel<nel/2+1;iel++)
    {
      if(parity==1) c.data[iel]=(data[iel]+data[(nel-iel)%nel])/2;
      else          c.data[iel]=(data[iel]-data[(nel-iel)%nel])/2;
    }

  return c;
}

TYPE constant_fit(VTYPE in,int tin,int tfin,const char *path=NULL)
{
  TYPE E(in.data[0]);

  E=0;
  double norm=0;
  for(int iel=max(tin,0);iel<=min(tfin,in.nel-1);iel++)
    {
      TYPE ele=in.data[iel];
      double err=in.data[iel].err();
      double weight=1/(err*err);
      if(!isnan(err))
	{
	  E+=ele*weight;
	  norm+=weight;
	}
    }
  E/=norm;
  
  if(path!=NULL)
    {
      double ym=E.med(),dy=E.err();
      ofstream out(path);
      out<<"@page size 800,600"<<endl;
      //error of the line
      out<<"@s0 line type 1"<<endl;      
      out<<"@s0 line color 7"<<endl;
      out<<"@s0 fill color 7"<<endl;
      out<<"@s0 fill type 1"<<endl;
      out<<tin<<" "<<ym-dy<<endl<<tfin<<" "<<ym-dy<<endl;
      out<<tfin<<" "<<ym+dy<<endl<<tin<<" "<<ym+dy<<endl;
      out<<tin<<" "<<ym-dy<<endl;
      out<<"&"<<endl;
      //central line
      out<<"@s1 line color 1"<<endl;
      out<<tin<<" "<<ym<<endl<<tfin<<" "<<ym<<endl;
      //plot the original data with error  
      out<<"&"<<endl;
      out<<"@type xydy"<<endl;      
      out<<"@s2 line type 0"<<endl;      
      out<<"@s2 symbol color 1"<<endl;
      out<<"@s2 errorbar color 1"<<endl;
      out<<"@s2 symbol 1"<<endl;
      out<<in;
      out.close();
    }
  
return E;
}

void linear_fit(VTYPE in,TYPE &m,TYPE &q,int tin,int tfin)
{
  double S,Sx,Sx2;
  TYPE Sxy(in.data[0]),Sy(in.data[0]);
  
  Sx2=S=Sx=0;
  Sxy=Sy=0;
  for(int iel=max(tin,0);iel<=min(tfin,in.nel-1);iel++)
    {
      double err=in.data[iel].err();
      double weight=1/(err*err);
      int x=iel;
      TYPE y=in.data[iel];

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

void linear_fit(TYPE &m,TYPE &q,double *x,VTYPE &y)
{
  double S,Sx,Sx2;
  TYPE Sxy(y.data[0]),Sy(y.data[0]);
  
  Sx2=S=Sx=0;
  Sxy=Sy=0;
  for(int iel=0;iel<y.nel;iel++)
    {
      double err=y.data[iel].err();
      double weight=1/(err*err);
      double xi=x[iel];
      TYPE yi=y.data[iel];
      
      S+=weight;
      Sx+=xi*weight;
      Sx2+=xi*xi*weight;
      Sxy+=xi*yi*weight;
      Sy+=yi*weight;
    }
  
  double delta=S*Sx2-Sx*Sx;
  m=(S*Sxy-Sx*Sy)/delta;
  q=(Sx2*Sy-Sxy*Sx)/delta;
}

VTYPE par_single_fun(double (*fun)(double,double*),VTYPE &x,VTYPE &par)
{
  int nx=x.nel;
  int npar=par.nel;
  
  VTYPE y(x);
  
  for(int ijack=0;ijack<x.njack+1;ijack++)
    {
      double dpar[npar];
      //create temporary double vector for pars
      for(int ipar=0;ipar<npar;ipar++) dpar[ipar]=par.data[ipar].data[ijack];
      for(int ix=0;ix<nx;ix++) y.data[ix].data[ijack]=fun(x.data[ix].data[ijack],dpar);
    }
  
  return y;
}
