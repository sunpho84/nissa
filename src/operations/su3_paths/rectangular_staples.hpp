#ifndef _RECTANGLE_STAPLE_HPP
#define _RECTANGLE_STAPLE_HPP

#include "squared_staples.hpp"

namespace nissa
{
  typedef squared_staples_t rectangular_staples_t;
  
  void compute_rectangular_staples_lx_conf(LxField<rectangular_staples_t>& out,
					   const LxField<quad_su3>& conf,
					   const LxField<squared_staples_t>& sq_staples);
  void compute_summed_rectangular_staples_lx_conf(LxField<quad_su3>& out,
						  const LxField<quad_su3>& conf,
						  const LxField<squared_staples_t>& sq_staples);
}

#endif
