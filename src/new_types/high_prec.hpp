#ifndef _HIGH_PREC_H
#define _HIGH_PREC_H

namespace nissa
{
  int high_prec_nbits();
  void float_high_prec_t_print(float_high_prec_t a);
  float_high_prec_t float_high_prec_t_pow_int(float_high_prec_t in,int d);
  float_high_prec_t float_high_prec_t_pow_int_frac(float_high_prec_t ext_in,int n,int d);
}

#endif
