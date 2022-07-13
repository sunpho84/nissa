#ifndef _ENDIANNESS_HPP
#define _ENDIANNESS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdint.h>

#ifndef EXTERN_ENDIANNESS
 #define EXTERN_ENDIANNESS extern
#endif

#include "routines/ios.hpp"

namespace nissa
{
  //types to revert uint16_t, float and double
  union uint16_t_reverter_t
  {
    uint16_t u;
    char c[2];
  };
  union float_reverter_t
  {
    float f;
    char c[4];
  };
  union double_reverter_t
  {
    double d;
    char c[8];
  };
  
  //endianness
  CUDA_MANAGED EXTERN_ENDIANNESS int little_endian;
  
  void check_endianness();
  
  namespace
  {
    template <typename T>
    CUDA_HOST_AND_DEVICE
    inline void swap(T& a,T& b)
    {
      T tmp=a;
      a=b;
      b=tmp;
    }
  }
  
  //revert the endianness of doubles
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(double *dest,double *sour,int ndoubles,int verbose)
  {
#ifndef COMPILING_FOR_DEVICE
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of %d doubles\n",ndoubles);
#endif
    
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
	double_reverter_t temp;
	temp.d=sour[idouble];
	swap(temp.c[0],temp.c[7]);
	swap(temp.c[1],temp.c[6]);
	swap(temp.c[2],temp.c[5]);
	swap(temp.c[3],temp.c[4]);
	dest[idouble]=temp.d;
      }
  }
  
  //revert the endianness of floats
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(float *dest,float *sour,int nfloats,int verbose)
  {
#ifndef COMPILING_FOR_DEVICE
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of %d floats\n",nfloats);
#endif
    
    for(int ifloat=0;ifloat<nfloats;ifloat++)
      {
	float_reverter_t temp;
	temp.f=sour[ifloat];
	swap(temp.c[0],temp.c[3]);
	swap(temp.c[1],temp.c[2]);
	dest[ifloat]=temp.f;
      }
  }
  
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(uint32_t *dest,uint32_t *sour,int nints,int verbose=1)
  {
    change_endianness((float*)dest,(float*)sour,nints,verbose);
  }
  
  CUDA_HOST_AND_DEVICE
  inline void change_endianness(uint16_t *dest,uint16_t *sour,int nshorts,int verbose=1)
  {
#ifndef COMPILING_FOR_DEVICE
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of %d uint16_t\n",nshorts);
#endif
    
    for(int ishort=0;ishort<nshorts;ishort++)
      {
	uint16_t_reverter_t temp;
	temp.u=sour[ishort];
	swap(temp.c[0],temp.c[1]);
	dest[ishort]=temp.u;
      }
  }
  
  ////////////////////Copy a vector of floats to doubles. Sweep is reversed to avoid overwriting////////////////
  
  //Do not change endianness
  inline void floats_to_doubles_same_endianness(double *dest,float *sour,int n,int verbose=1)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d floats to doubles\n",n);
    for(int i=n-1;i>=0;i--) dest[i]=(double)(sour[i]);
  }
  
  //Change endianness
  inline void floats_to_doubles_changing_endianness(double *dest,float *sour,int n,int verbose=1)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d floats to doubles changing endianness (two steps\n",n);
    floats_to_doubles_same_endianness(dest,sour,n,verbose);
    change_endianness(dest,dest,n,verbose);
  }
  
  ////////////////////Copy a vector of doubles to floats. Sweep is direct, to avoid overwriting////////////////
  
  //Do not change the endianness
  inline void doubles_to_floats_same_endianness(float *dest,double *sour,int n,int verbose=1)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d doubles to floats\n",n);
    for(int i=0;i<n;i++) dest[i]=(float)(sour[i]);
  }
  
  //Change endianness
  inline void doubles_to_floats_changing_endianness(float *dest,double *sour,int n,int verbose=1)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d doubles to floats changing endianness (two steps)\n",n);
    doubles_to_floats_same_endianness(dest,sour,n,verbose);
    change_endianness(dest,dest,n,verbose);
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
