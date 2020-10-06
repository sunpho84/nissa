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
  typedef mpf_class float_high_prec_t;
#elif HIGH_PREC == NATIVE_HIGH_PREC
  typedef float_256_class float_high_prec_t;
#else
 #error Unknwon high_prec: HIGH_PREC
#endif
  
  int high_prec_nbits();
  void float_high_prec_t_print(float_high_prec_t a);
  float_high_prec_t float_high_prec_t_pow_int(float_high_prec_t in,int d);
  float_high_prec_t float_high_prec_t_pow_int_frac(float_high_prec_t ext_in,int n,int d);
}

#endif
