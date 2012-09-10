#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/fourier.h"
#include "../routines/shift.h"

void mom_space_qq_vertex_function(spinspin v,int imom1,int imom2,quark_info qu,int mu)
{
  if(mu<0||mu>3) crash("mu=%d",mu);

  int summ=glblx_of_summ(imom1,imom2);

  double p=M_PI*(2*glb_coord_of_loclx[summ][mu]+qu.bc[mu])/glb_size[mu];
  double ph=p/2;

  spinspin_put_to_id(v);
  spinspin_prodassign_double(v,-sin(ph));
  spinspin_dirac_summ_the_prod_idouble(v,base_gamma+nissa_map_mu[mu],-cos(ph));
  spinspin_prodassign_double(v,glb_vol);
}
