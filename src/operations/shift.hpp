#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include <stdlib.h>
#include "new_types/su3.hpp"

namespace nissa
{
  void su3_vec_single_shift(su3 *u,int mu,int sign);
  void quad_su3_vec_single_shift(su3 *u,int mu,int sign);
  inline void su3_vec_single_shift(su3 *u,int signed_mu)
  {su3_vec_single_shift(u,((signed_mu<0)?-1:+1),labs(signed_mu));}
}

#endif
