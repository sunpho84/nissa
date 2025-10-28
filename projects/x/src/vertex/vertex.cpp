#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../propagators/twisted_propagator.hpp"
#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../routines/shift.hpp"

void mom_space_qq_vertex_function(spinspin v,double p1,double p2,quark_info qu,int mu)
{
  if(mu<0||mu>3) CRASH("mu=%d",mu);

  double s=(p1+p2)/2;
  spinspin_put_to_id(v);
  spinspin_prodassign_double(v,-sin(s)*glb_vol);
  spinspin_dirac_summ_the_prod_idouble(v,base_gamma+map_mu[mu],-cos(s)*glb_vol);
}
