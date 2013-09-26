#ifndef _ENDIANNESS_H
#define _ENDIANNESS_H

namespace nissa
{
  void check_endianness();
  void doubles_to_doubles_changing_endianness(double *dest,double *sour,int ndoubles,int verbose=1);
  void doubles_to_floats_changing_endianness(float *dest,double *sour,int n,int verbose=1);
  void doubles_to_floats_same_endianness(float *dest,double *sour,int n,int verbose=1);
  void floats_to_doubles_changing_endianness(double *dest,float *sour,int n,int verbose=1);
  void floats_to_doubles_same_endianness(double *dest,float *sour,int n,int verbose=1);
  void floats_to_floats_changing_endianness(float *dest,float *sour,int nfloats,int verbose=1);
  void uint16s_to_uint16s_changing_endianness(uint16_t *dest,uint16_t *sour,int nints,int verbose=1);
  void uint32s_to_uint32s_changing_endianness(uint32_t *dest,uint32_t *sour,int nints,int verbose=1);
  void uint64s_to_uint64s_changing_endianness(uint64_t *dest,uint64_t *sour,int nints,int verbose=1);
}

#endif
