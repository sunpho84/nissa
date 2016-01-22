#ifndef _BGQ_INTRINSIC_HPP
#define _BGQ_INTRINSIC_HPP

#if (defined BGQ) && (!defined BGQ_EMU)
  typedef vector4double reg_bi_complex;
#else
  typedef bi_complex reg_bi_complex;
#endif //BGQ_EMU

#include "intrinsic/declare.hpp"
#include "intrinsic/load.hpp"
#include "intrinsic/mergesplit.hpp"
#include "intrinsic/oper.hpp"
#include "intrinsic/store.hpp"

#endif
