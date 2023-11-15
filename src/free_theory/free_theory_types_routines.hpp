#ifndef _FREE_THEORY_TYPES_ROUTINES_HPP
#define _FREE_THEORY_TYPES_ROUTINES_HPP

#include "base/old_field.hpp"
#include "free_theory_types.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  gauge_info create_tlSym_gauge_info(const double& alpha,
				     const momentum_t& bc,
				     const double c1=-1.0/12);
  
  
  gauge_info create_Wilson_gauge_info(const double& alpha,
				      const momentum_t& bc);
  
  
  tm_quark_info create_twisted_quark_info(const double& kappa,
					  const double& mass,
					  const momentum_t& bc,
					  const int& r,
					  const double& zmp=0);
  
  
  tm_quark_info create_Wilson_quark_info(const double& kappa,
					 const momentum_t& bc);
  
  
  void get_spin_from_spinspin(LxField<spin0>& out,
			      const LxField<spinspin>& in,
			      const int& id_so);
  
  
  void put_spin_into_spinspin(LxField<spinspin>& out,
			      const LxField<spin0>& in,
			      const int& id_so);
}

#endif
