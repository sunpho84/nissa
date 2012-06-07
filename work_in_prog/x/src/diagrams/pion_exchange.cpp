#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../routines/fourier.h"

#include "propagator_self_energy.h"

/*
    __A__
   /  _} \
  0  {    X
   \__B__/
   
 */
 

void compute_twisted_propagator_exchange(spinspin *q_out,quark_info qu,gluon_info gl)
{
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}
