#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include <math.h>

namespace nissa
{
  double average_real_part_of_trace_of_rectangle_path(quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u);
  void average_trace_of_rectangle_path(complex tra,quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u);
  void su3_vec_single_shift(su3 *u,int mu,int sign);
  inline void su3_vec_single_shift(su3 *u,int signed_mu)
  {su3_vec_single_shift(u,((signed_mu<0)?-1:+1),std::abs(signed_mu));}
}

#endif
