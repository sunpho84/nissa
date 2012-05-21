#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/spin.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/random.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/operations/fft.h"

#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/fourier.h"
#include "../types/types.h"

void generate_stochastic_source_eta(spin1field *eta)
{
  //start (or restart) the random generator with the passed source
  if(!nissa_loc_rnd_gen_inited) crash("random gen not inited");
  
  //fill with Z4
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      comp_get_rnd(eta[ivol][mu],&(loc_rnd_gen[ivol]),RND_Z4);
  
  set_borders_invalid(eta);
}
