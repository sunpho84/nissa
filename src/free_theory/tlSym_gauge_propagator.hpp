#ifndef _TLSYM_GAUGE_PROPAGATOR_HPP
#define _TLSYM_GAUGE_PROPAGATOR_HPP

#include "free_theory_types.hpp"

namespace nissa
{
  void compute_mom_space_tlSym_gauge_propagator(spin1prop *prop,gauge_info gl);
  void compute_x_space_tlSym_gauge_propagator_by_fft(spin1prop *prop,gauge_info gl);
  void generate_stochastic_tlSym_gauge_propagator(spin1field *phi,spin1field *eta,gauge_info gl);
  void mom_space_tlSym_gauge_propagator_of_imom(spin1prop prop,gauge_info gl,int imom);
  void multiply_mom_space_tlSym_gauge_propagator(spin1field *out,spin1field *in,gauge_info gl);
  void multiply_x_space_tlSym_gauge_propagator_by_fft(spin1prop *out,spin1prop *in,gauge_info gl);
}

#endif
