#ifndef _FFT_HPP
#define _FFT_HPP

namespace nissa
{
  int bitrev(int in,int l2n);
  int find_max_pow2(int a);
  void data_coordinate_order_shift(complex *data,int ncpp,int mu0);
  void fft4d(complex *out,complex *in,int *dirs,int ncpp,double sign,int normalize);
  inline void fft4d(complex *out,complex *in,int ncpp,double sign,int normalize)
  {fft4d(out,in,all_dirs,ncpp,sign,normalize);}
}

#endif
