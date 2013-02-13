#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

double sqr(double a) {return a*a;}
int min_int(int a,int b) {if(a<b) return a;else return b;}
int max_int(int a,int b) {if(a>b) return a;else return b;}
double min_double(double a,double b) {if(a<b) return a;else return b;}
double max_double(double a,double b) {if(a>b) return a;else return b;}
void swap_doubles(double &d1,double &d2)
{
  double tmp=d1;
  d1=d2;
  d2=tmp;
}
void swap_ints(int &i1,int &i2)
{
  double tmp=i1;
  i1=i2;
  i2=tmp;
}

//return the log of the factorial of n
double lfact(double n)
{
  double log_factn=0;
  for(int i=1;i<=n;i++) log_factn+=log(i);
  return log_factn;
}

//return the log2 of N
int log2N(int N)
{
  int log2N=0;
  
  do log2N++;
  while ((2<<log2N)<N);
  
  return log2N;
}

//factorize a number
int factorize(int *list,int N)
{
  int nfatt=0;
  int fatt=2;

  while(N>1)
    {
      int div=N/fatt;
      int res=N-div*fatt;
      if(res!=0) fatt++;
      else 
	{
	  N=div;
	  list[nfatt]=fatt;
	  nfatt++;
	}
    }

  return nfatt;
}
