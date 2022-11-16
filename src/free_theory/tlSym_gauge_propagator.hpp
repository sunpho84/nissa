#ifndef _TLSYM_GAUGE_PROPAGATOR_HPP
#define _TLSYM_GAUGE_PROPAGATOR_HPP

#include "free_theory_types.hpp"
#include "new_types/spin.hpp"

namespace nissa
{
  bool zero_mode_subtraction_mask(const gauge_info& gl,const LocLxSite& imom);
  void compute_mom_space_tlSym_gauge_propagator(spin1prop *prop,const gauge_info& gl);
  void compute_x_space_tlSym_gauge_propagator_by_fft(spin1prop *prop,const gauge_info& gl);
  void generate_stochastic_tlSym_gauge_propagator_source(spin1field *eta);
  void generate_stochastic_tlSym_gauge_propagator(spin1field *phi,spin1field *eta,const gauge_info& gl);
  CUDA_HOST_AND_DEVICE void mom_space_tlSym_gauge_propagator_of_imom(spin1prop prop,const gauge_info& gl,int imom);
  void multiply_mom_space_tlSym_gauge_propagator(spin1field *out,spin1field *in,const gauge_info& gl);
  void multiply_x_space_tlSym_gauge_propagator_by_fft(spin1prop *out,spin1prop *in,const gauge_info& gl);
  void multiply_by_sqrt_tlSym_gauge_propagator(spin1field *out,spin1field *in,const gauge_info& gl);
  void multiply_by_tlSym_gauge_propagator(spin1field *out,spin1field *in,const gauge_info& gl);
  
  /// Compute the tadpole by taking the zero momentum ft of momentum prop
  Momentum compute_tadpole(const gauge_info& photon);
  double gluon_energy(gauge_info gl,const double virt,const GlbLxSite& imom);
}

#endif
