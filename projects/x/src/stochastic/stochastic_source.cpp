#include "../../../../src/new_types/new_types_definitions.hpp"
#include "../../../../src/new_types/spin.hpp"
#include "../../../../src/base/global_variables.hpp"
#include "../../../../src/base/random.hpp"
#include "../../../../src/base/debug.hpp"
#include "../../../../src/base/vectors.hpp"
#include "../../../../src/operations/fft.hpp"

#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../types/types.hpp"

void generate_stochastic_source_eta(spin1field *eta)
{
  //start (or restart) the random generator with the passed source
  if(!loc_rnd_gen_inited) CRASH("random gen not inited");
  
  //fill with Z4
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      comp_get_rnd(eta[ivol][mu],&(loc_rnd_gen[ivol]),RND_Z4);
  
  set_borders_invalid(eta);
}
