#ifndef _SHIFT_H
#define _SHIFT_H

namespace nissa
{
  double average_real_part_of_trace_of_rectangle_path(quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u);
  void average_trace_of_rectangle_path(complex tra,quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u);
  void su3_vec_single_shift(su3 *u,int mu,int sign);
}

#endif
