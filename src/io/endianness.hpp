#ifndef _ENDIANNESS_HPP
#define _ENDIANNESS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdint.h>

#ifndef EXTERN_ENDIANNESS
 #define EXTERN_ENDIANNESS extern
#endif

namespace nissa
{
  //endianness
  CUDA_MANAGED EXTERN_ENDIANNESS int little_endian;
  
  void check_endianness();
  void doubles_to_floats_changing_endianness(float *dest,double *sour,int n,int verbose=1);
  void doubles_to_floats_same_endianness(float *dest,double *sour,int n,int verbose=1);
  void floats_to_doubles_changing_endianness(double *dest,float *sour,int n,int verbose=1);
  void floats_to_doubles_same_endianness(double *dest,float *sour,int n,int verbose=1);
  
  CUDA_HOST_AND_DEVICE
  void change_endianness(float *dest,float *sour,int nfloats,int verbose=1);
  
  CUDA_HOST_AND_DEVICE
  void change_endianness(uint16_t *dest,uint16_t *sour,int nints,int verbose=1);
  
  CUDA_HOST_AND_DEVICE
  void change_endianness(double *dest,double *sour,int ndoubles,int verbose=1);
  
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(uint32_t *dest,uint32_t *sour,int nints,int verbose=1)
  {
    change_endianness((float*)dest,(float*)sour,nints,verbose);
  }
  
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(int *dest,int *sour,int nints,int verbose=1)
  {
    change_endianness((float*)dest,(float*)sour,nints,verbose);
  }
  
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(uint64_t *dest,uint64_t *sour,int nints,int verbose=1)
  {
    change_endianness((double*)dest,(double*)sour,nints,verbose);
  }
  
  template <class T>
  CUDA_HOST_AND_DEVICE
  void change_endianness(T &a,int verbose=0)
  {
    change_endianness(&a,&a,1,verbose);
  }
}

#undef EXTERN_ENDIANNESS

#endif
