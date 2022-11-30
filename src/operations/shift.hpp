#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include <stdlib.h>

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void su3_vec_single_shift(LxField<su3>& u,
			    const int& mu,
			    const int& sign);
}

#endif
