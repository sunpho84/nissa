#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "../base/global_variables.h"
#include "../base/random.h"

double sqr(double a) {return a*a;}
int min_int(int a,int b) {if(a<b) return a;else return b;}
int max_int(int a,int b) {if(a>b) return a;else return b;}
double min_double(double a,double b) {if(a<b) return a;else return b;}
double max_double(double a,double b) {if(a>b) return a;else return b;}

//swap two doubles
void swap_doubles(double &d1,double &d2)
{
  double tmp=d1;
  d1=d2;
  d2=tmp;
}

//swap two ints
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

//perform the metropolis test on passed value
int metro_test(double arg)
{
  double tresh=(arg<=0) ? 1 : exp(-arg);
  double toss=rnd_get_unif(&glb_rnd_gen,0,1);
  
  return toss<tresh;
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
