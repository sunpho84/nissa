#ifndef _BGQ_INTRINSIC_HPP
#define _BGQ_INTRINSIC_HPP

#include "new_types/complex.hpp"

namespace nissa
{
#if (defined BGQ) && (!defined BGQ_EMU)
  typedef vector4double reg_vir_complex;
#else
  typedef vir_complex reg_vir_complex;
#endif //BGQ_EMU
}

#include "intrinsic/declare.hpp"
#include "intrinsic/load.hpp"
#include "intrinsic/mergesplit.hpp"
#include "intrinsic/oper.hpp"
#include "intrinsic/store.hpp"

#endif
