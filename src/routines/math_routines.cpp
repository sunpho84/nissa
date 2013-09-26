#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/global_variables.hpp"
#include "base/random.hpp"

namespace nissa
{
  double sqr(double a) {return a*a;}
  
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
  
  //compute treshold
  double metro_tresh(double arg)
  {return (arg<=0) ? 1 : exp(-arg);}
  
  //perform the metropolis test on passed value
  int metro_test(double arg)
  {
    double tresh=metro_tresh(arg);
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
}
