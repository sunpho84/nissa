#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/random.hpp"
#include "new_types/complex.hpp"

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
  
  //recursive call - see below
  CUDA_HOST_AND_DEVICE void determinant(complex d,complex *m,int *s,int n,int N)
  {
    switch(n)
      {
      case 1:
	complex_copy(d,m[s[0]]);
	break;
      case 2:
	unsafe_complex_prod(  d,m[s[0]],m[N+s[1]]);
	complex_subt_the_prod(d,m[s[1]],m[N+s[0]]);
	break;
      default:
	complex_put_to_zero(d);
	for(int p=0;p<n;p++)
	  {
	    //prepare submatrix
	    int *l=new int[n-1];
	    for(int i=0;i<p;i++) l[i]=s[i];
	    for(int i=p;i<n-1;i++) l[i]=s[i+1];
	    
	    //compute determinant
	    complex in_det;
	    determinant(in_det,m+N,l,n-1,N);
	    
	    //summ or subtract the product
	    void (*fun[2])(complex,const complex,const complex)={complex_summ_the_prod,complex_subt_the_prod};
	    fun[p%2](d,m[s[p]],in_det);
	    
	    delete[] l;
	  }
      }
  }
  
  //compute the determinant of a NxN matrix through a recursive formula or eigen
  CUDA_HOST_AND_DEVICE void matrix_determinant(complex d,complex *m,int n)
  {
#if USE_EIGEN
    using cpp_complex=std::complex<double>;
    using cpp_matrix=cpp_complex*;
    using eig_matrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using map_eig_matrix=Eigen::Map<eig_matrix>;
    
    cpp_complex& cpp_d=*reinterpret_cast<cpp_complex*>(d);
    cpp_matrix cpp_m=reinterpret_cast<cpp_complex*>(m);
    map_eig_matrix eig_m(cpp_m,n,n);
    
    cpp_d=eig_m.determinant();
#else
    int *l=new int[n];
    for(int i=0;i<n;i++) l[i]=i;
    determinant(d,m,l,n,n);
    delete[] l;
#endif
  }
}
