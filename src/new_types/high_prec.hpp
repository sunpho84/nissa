#ifndef _HIGH_PREC_HPP
#define _HIGH_PREC_HPP

#if HIGH_PREC == GMP_HIGH_PREC
 #include <gmpxx.h>
#endif

#if HIGH_PREC == NATIVE
 #include "new_types/float_256.hpp"
#endif

namespace nissa
{
#if HIGH_PREC == GMP_HIGH_PREC
  
  //uses mpf for high precision
  typedef mpf_class float_high_prec_t;
  
  //stores the number of bit oh mpf precision
  extern int mpf_precision;
  
  //default number of bits
  #define NISSA_DEFAULT_MPF_PRECISION 256
  
#elif HIGH_PREC == NATIVE_HIGH_PREC
  
  //revert to native implementation
  typedef float_256_class float_high_prec_t;
  
#else
 #error Unknwon high_prec: HIGH_PREC
#endif
  
  void init_high_precision();
  int high_prec_nbits();
  void float_high_prec_t_print(float_high_prec_t a);
  float_high_prec_t float_high_prec_t_pow_int(float_high_prec_t in,int d);
  float_high_prec_t float_high_prec_t_pow_int_frac(float_high_prec_t ext_in,int n,int d);
}

#endif
