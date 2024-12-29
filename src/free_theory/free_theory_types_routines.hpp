#ifndef _FREE_THEORY_TYPES_ROUTINES_HPP
#define _FREE_THEORY_TYPES_ROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
#include <free_theory/free_theory_types.hpp>

namespace nissa
{
  gauge_info create_tlSym_gauge_info(const gauge_info::which_gauge_t& which_gauge,
				     const Momentum& bc,
				     const double c1=-1.0/12);
  
  gauge_info create_Wilson_gauge_info(const gauge_info::which_gauge_t& which_gauge,
				      const Momentum& bc);
  
  TmQuarkInfo create_twisted_quark_info(const double& kappa,
					const double& mass,
					const Momentum& bc,
					const int& r,
					const double& zmp=0);
  
  TmQuarkInfo create_Wilson_quark_info(const double& kappa,
				       const Momentum& bc);
  
  void get_spin_from_spinspin(LxField<spin>& out,
			      const LxField<spinspin>& in,
			      const int& id_so);
  
  void put_spin_into_spinspin(LxField<spinspin>& out,
			      const LxField<spin>& in,
			      const int& id_so);
}

#endif
