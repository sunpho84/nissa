#ifndef _TLSYM_GAUGE_PROPAGATOR_HPP
#define _TLSYM_GAUGE_PROPAGATOR_HPP

#include "base/field.hpp"
#include "free_theory_types.hpp"
#include "new_types/spin.hpp"

namespace nissa
{
  bool zero_mode_subtraction_mask(const gauge_info& gl,
				  const int& imom);
  
  void compute_mom_space_tlSym_gauge_propagator(LxField<spin1prop>& prop,
						const gauge_info& gl);
  
  void compute_x_space_tlSym_gauge_propagator_by_fft(LxField<spin1prop>& prop,
						     const gauge_info& gl);
  
  void generate_stochastic_tlSym_gauge_propagator_source(LxField<spin1field>& eta);
  
  void generate_stochastic_tlSym_gauge_propagator(LxField<spin1field>& phi,
						  LxField<spin1field>& eta,
						  const gauge_info& gl);
  
  CUDA_HOST_AND_DEVICE void mom_space_tlSym_gauge_propagator_of_imom(spin1prop& prop,
								     const gauge_info& gl,
								     const int& imom);
  
  void multiply_mom_space_tlSym_gauge_propagator(LxField<spin1field>& out,
						 const LxField<spin1field>& in,
						 const gauge_info& gl);
  
  void multiply_x_space_tlSym_gauge_propagator_by_fft(LxField<spin1prop>& out,
						      const LxField<spin1prop>& in,
						      const gauge_info& gl);
  
  void multiply_by_sqrt_tlSym_gauge_propagator(LxField<spin1field>& out,
					       const LxField<spin1field>& in,
					       const gauge_info& gl);
  
  void multiply_by_tlSym_gauge_propagator(LxField<spin1field>& out,
					  const LxField<spin1field>& in,
					  const gauge_info& gl);
  
  Momentum compute_tadpole(const gauge_info& photon);
  
  double gluon_energy(const gauge_info& gl,
		      const double& virt,
		      const int& imom);
}

#endif
