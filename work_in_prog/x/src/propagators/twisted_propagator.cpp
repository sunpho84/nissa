#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/geometry/geometry_lx.h"
#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/complex.h"
#include "../../../../src/new_types/spin.h"
#include "../../../../src/operations/fft.h"
#include "../../../../src/operations/gaugeconf.h"
#include "../../../../src/inverters/twisted_mass/cg_invert_tmDeoimpr.h"

#include "../types/types.h"
#include "../routines/fourier.h"
#include "../inverters/cg_eoprec_twisted_free_operator.h"

void mom_space_twisted_propagator_of_imom(spin1prop prop,quark_info qu,int imom)
{
  double kappa=qu.kappa;
  double mass=qu.mass;
  
  double M=1/(2*kappa)-4;
  
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
      complex_prod_double(          prop[ig][base_gamma[0].pos[ig]],base_gamma[0].entr[ig],Mp*rep_den);
      complex_summ_the_prod_idouble(prop[ig][base_gamma[1].pos[ig]],base_gamma[1].entr[ig],-sin_mom[1]*rep_den);
      complex_summ_the_prod_idouble(prop[ig][base_gamma[2].pos[ig]],base_gamma[2].entr[ig],-sin_mom[2]*rep_den);
      complex_summ_the_prod_idouble(prop[ig][base_gamma[3].pos[ig]],base_gamma[3].entr[ig],-sin_mom[3]*rep_den);
      complex_summ_the_prod_idouble(prop[ig][base_gamma[4].pos[ig]],base_gamma[4].entr[ig],-sin_mom[0]*rep_den);
      complex_summ_the_prod_idouble(prop[ig][base_gamma[5].pos[ig]],base_gamma[5].entr[ig],-mass*rep_den);
    }
}

//compute the twisted quark propagator in the momentum space
void multiply_mom_space_twisted_propagator(spin *out,spin *in,quark_info qu)
{
  nissa_loc_vol_loop(imom)
    {
      spin1prop prop;
      mom_space_twisted_propagator_of_imom(prop,qu,imom);
      safe_spinspin_spin_prod(out[imom],prop,in[imom]);
    }
  
  set_borders_invalid(out);
}

//compute the twisted quark propagator in the momentum space
void compute_mom_space_twisted_propagator(spinspin *prop,quark_info qu)
{
  nissa_loc_vol_loop(imom)
    mom_space_twisted_propagator_of_imom(prop[imom],qu,imom);
  
  set_borders_invalid(prop);
}

//pass from p to x space
void compute_x_space_twisted_propagator_by_fft(spinspin *prop,quark_info qu)
{
  compute_mom_space_twisted_propagator(prop,qu);
  pass_spinspin_from_mom_to_x_space(prop,prop,qu.bc);
}

//multiply the source for the twisted propagator by inverting twisted Dirac operator
void multiply_x_space_twisted_propagator_by_inverting(spin *prop,spin *ext_source,quark_info qu)
{inv_tmD_cg_eoprec_eos(prop,NULL,qu,1000000,1.e-28,ext_source);}
void multiply_x_space_twisted_propagator_by_inverting(spinspin *prop,spinspin *ext_source,quark_info qu)
{
  //source and temp prop
  spin *tsource=nissa_malloc("tsource",loc_vol+bord_vol,spin);
  spin *tprop=nissa_malloc("tprop",loc_vol,spin);
    
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      //prepare the source
      nissa_loc_vol_loop(ivol)
	for(int id_si=0;id_si<4;id_si++)
	  memcpy(tsource[ivol][id_si],ext_source[ivol][id_si][id_so],sizeof(complex));
      
      //invert and copy into the spinspin
      multiply_x_space_twisted_propagator_by_inverting(tprop,tsource,qu);
      nissa_loc_vol_loop(ivol)
	for(int id_si=0;id_si<4;id_si++)
	  memcpy(prop[ivol][id_si][id_so],tprop[ivol][id_si],sizeof(complex));
    }
  
  set_borders_invalid(prop);
  
  nissa_free(tsource);
  nissa_free(tprop);
}

void compute_x_space_twisted_propagator_by_inverting(spinspin *prop,quark_info qu)
{
  spinspin *delta=nissa_malloc("delta",loc_vol+bord_vol,spinspin);
  memset(delta,0,sizeof(spinspin)*loc_vol);
  if(rank==0) for(int id=0;id<4;id++) delta[0][id][id][0]=1;
  
  multiply_x_space_twisted_propagator_by_inverting(prop,delta,qu);
  
  set_borders_invalid(prop);
  
  nissa_free(delta);
}
