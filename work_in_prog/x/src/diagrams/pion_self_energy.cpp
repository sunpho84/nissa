#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../routines/fourier.h"

#include "propagator_self_energy.h"

void compute_twisted_propagator_with_self_energy_insertion(spinspin *q_out,quark_info qu,gluon_info gl)
{
  compute_self_energy_twisted_propagator_in_x_space(q_out,qu,gl);
  pass_spinspin_from_x_to_mom_space(q_out,q_out,qu.bc);
  
  nissa_loc_vol_loop(imom)
    {
      spinspin s,t;
      mom_space_twisted_propagator_of_imom(s,qu,imom);
      unsafe_spinspin_spinspin_prod(t,s,q_out[imom]);
      unsafe_spinspin_spinspin_prod(q_out[imom],s,t);
    }
  
  pass_spinspin_from_mom_to_x_space(q_out,q_out,qu.bc);
}
