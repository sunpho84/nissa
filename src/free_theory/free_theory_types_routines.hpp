#ifndef _FREE_THEORY_TYPES_ROUTINES_HPP
#define _FREE_THEORY_TYPES_ROUTINES_HPP

#include "new_types/su3.hpp"
#include "free_theory_types.hpp"

namespace nissa
{
  gauge_info create_tlSym_gauge_info(const double& alpha,const Momentum& bc,const double& c1=-1.0/12);
  gauge_info create_Wilson_gauge_info(const double& alpha,const Momentum& bc);
  tm_quark_info create_Wilson_quark_info(const double& kappa,const Momentum& bc);
  tm_quark_info create_twisted_quark_info(const double& kappa,const double& mass,const Momentum& bc,int r,const double& zmp=0);
  void get_spin_from_spinspin(spin *out,spinspin *in,const int& id_so);
  void put_spin_into_spinspin(spinspin *out,spin *in,const int& id_so);
}

#endif
