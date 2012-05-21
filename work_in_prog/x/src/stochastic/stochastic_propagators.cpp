#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/spin.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/operations/fft.h"

#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/fourier.h"
#include "../types/types.h"

#include "stochastic_source.h"

//perform the multiplication in the momentum space in order to avoid convolution
void generate_stochastic_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl)
{
  //pass to mom space
  pass_spin1field_from_x_to_mom_space(phi,eta,gl.bc);
  
  //multiply by prop
  multiply_mom_space_tlSym_gluon_propagator(phi,phi,gl);
  
  //takes the anti-fast fourier transform of eta
  pass_spin1field_from_mom_to_x_space(phi,phi,gl.bc);
  
  set_borders_invalid(phi);
}

void generate_stochastic_source_and_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl)
{
  generate_stochastic_source_eta(eta);
  generate_stochastic_tlSym_gluon_propagator(phi,eta,gl);
}
