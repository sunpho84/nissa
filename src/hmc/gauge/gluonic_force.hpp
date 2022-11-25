#ifndef _GLUONIC_FORCE_HPP
#define _GLUONIC_FORCE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "hmc/theory_pars.hpp"

namespace nissa
{
  void gluonic_force_finish_computation(LxField<quad_su3>& F,
					const LxField<quad_su3>& conf);
  
  void compute_gluonic_force_lx_conf_do_not_finish(LxField<quad_su3>& F,
						   const LxField<quad_su3>& conf,
						   const theory_pars_t& physics);
  
  void compute_gluonic_force_lx_conf(LxField<quad_su3>& F,
				     const LxField<quad_su3>& conf,
				     const theory_pars_t& physics);
}

#endif
