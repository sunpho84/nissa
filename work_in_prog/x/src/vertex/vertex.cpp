#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/fourier.h"
#include "../routines/shift.h"

void mom_space_qq_vertex_function(spinspin v,double p1,double p2,quark_info qu,int mu)
{
  if(mu<0||mu>3) crash("mu=%d",mu);

  double s=(p1+p2)/2;
  spinspin_put_to_id(v);
  spinspin_prodassign_double(v,-sin(s)*glb_vol);
  spinspin_dirac_summ_the_prod_idouble(v,base_gamma+nissa_map_mu[mu],-cos(s)*glb_vol);
}
