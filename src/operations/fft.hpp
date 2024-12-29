#ifndef _FFT_HPP
#define _FFT_HPP

#include "base/field.hpp"
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
	     const WhichDirs& dirs,
	     const int& ncpp,
	     const double& sign,
	     const bool& normalize);
  
  inline void fft4d(complex *out,
		     const complex *in,
		     const int& ncpp,
		     const double& sign,
		     const bool& normalize)
   {
     fft4d(out,in,allDirs,ncpp,sign,normalize);
   }
  
  inline void fft4d(void *out,
		    const void* in,
		    const int& ncpp,
		    const double& sign,
		    int normalize)
  {
    fft4d((complex*)out,(const complex*)in,allDirs,ncpp,sign,normalize);
  }
  
  // template <class T>
  // void fft4d(T *x,double sign,int normalize)
  // {fft4d(x,x,sign,normalize);}
  
  inline void fft4d(void *x,
		    const int& ncpp,
		    const double& sign,
		    const bool normalize)
  {
    fft4d(x,x,ncpp,sign,normalize);
  }
  
  /// Fourier transform of the field
  template <typename T,
	    FieldLayout FL>
  void fft4d(LxField<T,FL>& f,
	     const double& sign,
	     const bool& normalize)
  {
    static_assert(LxField<T>::nInternalDegs%2==0,"Needs to be able to interpret the field as complex");
    
    if constexpr(FL!=FieldLayout::CPU)
      {
	LxField<T,FieldLayout::CPU> tmp("tmp");
	tmp=f;
	master_printf("Making a temporary copy of the field to be 4D-fourier-transformed\n");
	
	fft4d(tmp,sign,normalize);
	
	f=tmp;
      }
    else
      fft4d(f._data,LxField<T>::nInternalDegs/2,sign,normalize);
  }
}

#endif
