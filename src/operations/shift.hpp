#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include <stdlib.h>

#include <geometry/geometry_lx.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void su3_vec_single_shift(su3 *u,const Direction& mu,int sign);
}

#endif
