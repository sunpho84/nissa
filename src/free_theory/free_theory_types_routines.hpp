#ifndef _FREE_THEORY_TYPES_ROUTINES_HPP
#define _FREE_THEORY_TYPES_ROUTINES_HPP

#include "new_types/su3.hpp"
#include "free_theory_types.hpp"

namespace nissa
{
  gauge_info create_tlSym_gauge_info(double alpha,const momentum_t& bc,double c1=-1.0/12);
  gauge_info create_Wilson_gauge_info(double alpha,const momentum_t& bc);
  tm_quark_info create_Wilson_quark_info(double kappa,const momentum_t& bc);
  tm_quark_info create_twisted_quark_info(double kappa,double mass,const momentum_t& bc,int r,double zmp=0);
  void get_spin_from_spinspin(spin *out,spinspin *in,int id_so);
  void put_spin_into_spinspin(spinspin *out,spin *in,int id_so);
}

#endif
