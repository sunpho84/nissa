#ifndef _FFT_HPP
#define _FFT_HPP

#include "base/old_field.hpp"
// #include "new_types/complex.hpp"

namespace nissa
{
  enum{FFT_NO_NORMALIZE=0,FFT_NORMALIZE=1};
#define FFT_PLUS +1.0
#define FFT_MINUS -1.0
  
  // int bitrev(int in,int l2n);
  // int find_max_pow2(int a);
  // void data_coordinate_order_shift(complex *data,int ncpp,int mu0);
  
  void fft4d(complex *out,
	     const complex *in,
	     const which_dir_t& dirs,
	     const int& ncpp,
	     const double& sign,
	     const bool& normalize);
  
  inline void fft4d(complex *out,
		     const complex *in,
		     const int& ncpp,
		     const double& sign,
		     const bool& normalize)
   {
     fft4d(out,in,all_dirs,ncpp,sign,normalize);
   }
  
  template <typename T>
  void fft4d(T *out,
	     const T *in,
	     const double& sign,
	     int normalize)
  {
    fft4d((complex*)out,(const complex*)in,all_dirs,sizeof(T)/sizeof(complex),sign,normalize);
  }
  
  // template <class T>
  // void fft4d(T *x,double sign,int normalize)
  // {fft4d(x,x,sign,normalize);}
  
  template <typename T>
  void fft4d(T *x,
	     const double& sign,
	     const bool normalize)
  {
    fft4d(x,x,sign,normalize);
  }
  
  /// Fourier transform of the field
  template <typename T,
	    FieldLayout FL>
  void fft4d(LxField<T,FL>& f,
	const double& sign,
	const bool& normalize)
  {
    LxField<T,FieldLayout::CPU> tmp("tmp");
    tmp=f;
    
    fft4d(tmp,sign,normalize);
    
    f=tmp;
  }
  
  template <typename T>
  void fft4d(LxField<T,FieldLayout::CPU>& f,
	const double& sign,
	const bool& normalize)
  {
    fft4d(f._data,sign,normalize);
  }
}

#endif
