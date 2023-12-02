#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#if USE_EIGEN
 #include <unsupported/Eigen/MatrixFunctions>
#endif

namespace nissa
{
  //return the log of the factorial of n
  double lfact(double n)
  {
    double log_factn=0;
    for(int i=1;i<=n;i++) log_factn+=log(i);
    return log_factn;
  }
  
  //return the bit inverse of an int
  int bitrev(int in,int l2n)
  {
    int out=0;
    
    for(int i=0;i<l2n;i++) if(in & (1<<i)) out+=(1<<(l2n-i-1));
    
    return out;
  }
  
  //return the powers of n contained in the input
  int find_max_pow2(int a)
  {
    int nl=0;
    while((a&0x1)==0)
      {
	nl++;
	a>>=1;
      }
    
    return nl;
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
