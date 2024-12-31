#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.hpp"
#include "../../../../src/base/debug.hpp"
#include "../../../../src/new_types/new_types_definitions.hpp"
#include "../../../../src/new_types/spin.hpp"
#include "../../../../src/new_types/dirac.hpp"

#include "../types/types.hpp"

void mom_space_twisted_propagator_g2_d2_corr_of_imom(spin1prop prop,quark_info qu,gluon_info gl,int imom)
{
  if(qu.kappa!=1.0/8 || qu.mass!=0) CRASH("implemented only for massless quark");
  if(gl.c1!=0 && gl.c1!=-1.0/12) CRASH("implemented only for Wilson or Symanzik (c1=-1/12) gluons");
  double lambda=gl.alpha;
  
  //implement equation 26 of 0907.0381
  
  //Table VII
  double et2_W[3]={7.4696262,8.2395316,-3.21623339};
  double et2_I[3]={5.95209802,7.16084231,-3.0693492};
  double *et2;
  if(gl.c1==0) et2=et2_W;
  else         et2=et2_I;
  
  momentum_t p,p3;
  double p2=0;
  for(int mu=0;mu<4;mu++)
    {
      p[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
      p3[mu]=p[mu]*p[mu]*p[mu];
      p2+=p[mu]*p[mu];
    }
  
  double c_id=et2[0]-5.39301570*lambda-0.5*(3-2*lambda)*log(p2);
  
  memset(prop,0,sizeof(spin1prop));
  
  spinspin_dirac_summ_the_prod_double(prop,&(base_gamma[0]),c_id);
}
