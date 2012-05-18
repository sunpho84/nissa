#include <string.h>
#include <math.h>

#include "../../src/base/global_variables.h"
#include "../../src/base/vectors.h"
#include "../../src/geometry/geometry_lx.h"
#include "../../src/new_types/new_types_definitions.h"
#include "../../src/new_types/complex.h"
#include "../../src/new_types/su3.h"
#include "../../src/operations/fft.h"
#include "../../src/operations/gaugeconf.h"
#include "../../src/inverters/twisted_mass/cg_invert_tmDeoimpr.h"

#include "types.h"
#include "fourier.h"

//compute the twisted quark propagator in the momentum space
void compute_mom_space_twisted_propagator(spinspin *prop,quark_info qu)
{
  double kappa=qu.kappa;
  double mass=qu.mass;
  
  double M=1/(2*kappa)-4;
  
  //reset the propagator
  memset(prop,0,loc_vol*sizeof(spinspin));
  
  nissa_loc_vol_loop(imom)
  {
    momentum_t sin_mom;
    double Mp=0,sin2_mom=0;
    for(int mu=0;mu<4;mu++)
      {
	double p=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
	double ph=p/2;
	sin_mom[mu]=sin(p);
	double sinph=sin(ph);
	Mp+=sinph*sinph;
	sin2_mom+=sin_mom[mu]*sin_mom[mu];
      }
    Mp=2*Mp+M;
      
    double rep_den=1/(sin2_mom+Mp*Mp+mass*mass);
      
    for(int ig=0;ig<4;ig++)
      {
	complex_prod_double(          prop[imom][ig][base_gamma[0].pos[ig]],base_gamma[0].entr[ig],Mp*rep_den);
	complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[1].pos[ig]],base_gamma[1].entr[ig],-sin_mom[1]*rep_den);
	complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[2].pos[ig]],base_gamma[2].entr[ig],-sin_mom[2]*rep_den);
	complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[3].pos[ig]],base_gamma[3].entr[ig],-sin_mom[3]*rep_den);
	complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[4].pos[ig]],base_gamma[4].entr[ig],-sin_mom[0]*rep_den);
	complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[5].pos[ig]],base_gamma[5].entr[ig],-mass*rep_den);
      }
  }
}

//pass from p to x space
void compute_x_space_twisted_propagator_by_fft(spinspin *prop,quark_info qu)
{
  compute_mom_space_twisted_propagator(prop,qu);
  pass_spinspin_from_mom_to_x_space(prop,prop,qu.bc);
}

//compute numerically the twisted propagator by making use of the already implmented inverters
void compute_x_space_twisted_propagator_by_inverting(spinspin *prop,quark_info qu)
{
  //source and temp prop
  spincolor *source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  spincolor *tprop=nissa_malloc("tprop",loc_vol,spincolor);
  
  //identical configuration+boundaries
  quad_su3 *conf=nissa_malloc("idconf",loc_vol+bord_vol,quad_su3);
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      su3_put_to_id(conf[ivol][mu]);
  
  //put boundary condition
  put_boundaries_conditions(conf,qu.bc,0,0);
  
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      //prepare the source
      memset(source,0,sizeof(spincolor)*loc_vol);  
      if(rank==0) source[0][id_so][0][0]=1;
      
      //invert and copy into the spinspin
      inv_tmD_cg_eoprec_eos(tprop,NULL,conf,qu.kappa,qu.mass,1000000,1.e-28,source);
      nissa_loc_vol_loop(ivol)
	for(int id_si=0;id_si<4;id_si++)
	  memcpy(prop[ivol][id_si][id_so],tprop[ivol][id_si][0],sizeof(complex));
    }
  
  //remove the boundary condition
  nissa_loc_vol_loop(ivol)
    {
      //compute the phase
      double arg=0;
      for(int mu=0;mu<4;mu++) arg+=qu.bc[mu]*glb_coord_of_loclx[ivol][mu]/glb_size[mu];
      arg*=M_PI;
      complex phase={cos(arg),sin(arg)};
      
      //multiply for the phase
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(prop[ivol][id1][id2],prop[ivol][id1][id2],phase);
    }
  
  nissa_free(conf);
  nissa_free(source);
  nissa_free(tprop);
}
