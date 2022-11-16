#ifndef _FFT_HPP
#define _FFT_HPP

#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"

namespace nissa
{
  enum{FFT_NO_NORMALIZE=0,FFT_NORMALIZE=1};
#define FFT_PLUS  +1.0
#define FFT_MINUS -1.0
  
  int bitrev(int in,int l2n);
  int find_max_pow2(int a);
  void data_coordinate_order_shift(complex *data,int ncpp,int mu0);
  void fft4d(complex *out,complex *in,const Coords<bool>& dirs,int ncpp,double sign,int normalize);
  inline void fft4d(complex *out,complex *in,int ncpp,double sign,int normalize)
  {
    fft4d(out,in,all_dirs,ncpp,sign,normalize);
  }
  
  template <class T>
  void fft4d(T *out,T *in,double sign,int normalize)
  {
    fft4d((complex*)out,(complex*)in,all_dirs,sizeof(T)/sizeof(complex),sign,normalize);
  }
  
  template <class T>
  void fft4d(T *x,double sign,int normalize)
  {
    fft4d(x,x,sign,normalize);
  }
}

#endif
