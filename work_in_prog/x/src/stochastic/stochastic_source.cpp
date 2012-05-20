#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/spin.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/random.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/operations/fft.h"

#include "../propagators/tlSym_gluon_propagator.h"
#include "../types/types.h"

void generate_stochastic_source_eta(spin1field *eta,int seed)
{
  //start (or restart) the random generator with the passed source
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  start_glb_rnd_gen(seed);
  
  //fill with Z4
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      comp_get_rnd(eta[ivol][mu],&(loc_rnd_gen[ivol]),RND_Z4);
  
  set_borders_invalid(eta);
}

//perform the multiplication in the momentum space in order to avoid convolution
void generate_stochastic_tlSym_gluon_propagator(spin1field *phi,spin1field *eta,gluon_info gl)
{
  crash("bc still need to be implemented");
  
  //takes the fast fourier transform of eta
  fft4d((complex*)phi,(complex*)eta,4,-1,0);
  
  //computes the propagator in mom space
  spin1prop *gl_mom_prop=nissa_malloc("gl_mom_prop",loc_vol,spin1prop);
  compute_mom_space_tlSym_gluon_propagator(gl_mom_prop,gl);
  
  //product (spin1prop is an alias for spinspin)
  nissa_loc_vol_loop(imom)
    safe_spinspin_spin_prod((phi[imom]),gl_mom_prop[imom],phi[imom]);
  
  //takes the anti-fast fourier transform of eta
  fft4d((complex*)phi,(complex*)phi,4,1,1);
  
  nissa_free(gl_mom_prop);
}
