#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include <stdlib.h>
#include "new_types/su3.hpp"

namespace nissa
{
  void su3_vec_single_shift(su3 *u,int mu,int sign);
  void quad_su3_vec_single_shift(su3 *u,int mu,int sign);
}

#endif
