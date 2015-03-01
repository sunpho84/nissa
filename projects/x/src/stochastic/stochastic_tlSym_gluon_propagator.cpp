#include "nissa.hpp" 
using namespace nissa;

#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../types/types.hpp"

#include "stochastic_source.hpp"

//perform the multiplication in the momentum space in order to avoid convolution
void generate_stochastic_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl)
{
  //pass to mom space
  pass_spin1field_from_x_to_mom_space(phi,eta,gl.bc);
  
  //multiply by prop
  multiply_mom_space_tlSym_gluon_propagator(phi,phi,gl);
  
  //put volume normalization due to convolution
  NISSA_LOC_VOL_LOOP(imom)
    spin_prodassign_double(phi[imom],glb_vol);  
  
  //takes the anti-fast fourier transform of eta
  pass_spin1field_from_mom_to_x_space(phi,phi,gl.bc);
  
  set_borders_invalid(phi);
}

void generate_stochastic_source_and_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl)
{
  generate_stochastic_source_eta(eta);
  generate_stochastic_tlSym_gluon_propagator(phi,eta,gl);
}
