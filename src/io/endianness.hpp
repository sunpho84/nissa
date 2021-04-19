#ifndef _ENDIANNESS_HPP
#define _ENDIANNESS_HPP

#include <stdint.h>

#ifndef EXTERN_ENDIANNESS
 #define EXTERN_ENDIANNESS extern
#endif

namespace nissa
{
  //endianness
  EXTERN_ENDIANNESS int little_endian;
  
  void check_endianness();
  void doubles_to_floats_changing_endianness(float* dest,const double* sour,const int& n,const bool& verbose=1);
  void doubles_to_floats_same_endianness(float* dest,const double* sour,const int& n,const bool& verbose=1);
  void floats_to_doubles_changing_endianness(double* dest,const float* sour,const int& n,const bool& verbose=1);
  void floats_to_doubles_same_endianness(double* dest,const float* sour,const int& n,const bool& verbose=1);
  void change_endianness(float* dest,const float* sour,const int& nfloats,const bool& verbose=1);
  void change_endianness(uint16_t* dest,const uint16_t* sour,const int& nints,const bool& verbose=1);
  void change_endianness(double* dest,const double* sour,const int& ndoubles,const bool& verbose=1);
  
  inline void change_endianness(uint32_t* dest,const uint32_t* sour,const int& nints,const bool& verbose=1)
  {
    change_endianness((float*)dest,(const float*)sour,nints,verbose);
  }
  
  inline void change_endianness(int* dest,const int* sour,const int& nints,const bool& verbose=1)
  {
    change_endianness((float*)dest,(const float*)sour,nints,verbose);
  }
  
  inline void change_endianness(uint64_t* dest,const uint64_t* sour,const int& nints,const bool& verbose=1)
  {
    change_endianness((double*)dest,(const double*)sour,nints,verbose);
  }
  
  inline void change_endianness(int64_t* dest,const int64_t* sour,const int& nints,const bool& verbose=1)
  {
    change_endianness((double*)dest,(const double*)sour,nints,verbose);
  }
  
  template <class T> void change_endianness(T &a,const bool& verbose=0)
  {
    change_endianness(&a,&a,1,verbose);
  }
}

#undef EXTERN_ENDIANNESS

#endif
