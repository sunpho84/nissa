#include <string.h>
#include <math.h>

#include "../../../../src/nissa.hpp"
using namespace std;

#include "../types/types.hpp"
#include "../inverters/cg_Wilson_gluon_operator.hpp"

//compute the Wilson action gluon propagator in the momentum space according to P.Weisz
void compute_mom_space_Wilson_gluon_propagator(spin1prop *prop,gluon_info gl)
{
  //check absence of zero modes
  int zmpres=1;
  for(int mu=0;mu<4;mu++) zmpres&=(gl.bc[mu]==0);
  //if(zmpres) crash("zero mode present, prop not defined");
  
  //reset the propagator
  memset(prop,0,loc_vol*sizeof(spin1prop));
  
  NISSA_LOC_VOL_LOOP(imom)
    {
      //momentum
      momentum_t k,kt;
      double kt2=0;
      for(int mu=0;mu<4;mu++)
	{
	  k[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+gl.bc[mu])/glb_size[mu]; //lattice momentum
	  kt[mu]=2*sin(k[mu]/2);
	  kt2+=kt[mu]*kt[mu];
	}
      
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    if(kt2!=0)
	      {
		if(mu==nu) prop[imom][mu][nu][RE]=1;
		prop[imom][mu][nu][RE]-=(1-gl.alpha)*kt[mu]*kt[nu]/kt2;
		prop[imom][mu][nu][RE]/=kt2;
	      }
	    else
	      prop[imom][mu][nu][RE]=prop[imom][mu][nu][IM]=0;
      
	    complex_prodassign_double(prop[imom][mu][nu],1.0/glb_vol);
	  }
    }
}

//compute the Wilson action gluon propagator in the x space by taking the fft of that in momentum space
void compute_x_space_Wilson_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl)
{
  compute_mom_space_Wilson_gluon_propagator(prop,gl);
  pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc);
}

//multiply by the Wilson action gluon propagator in the x space by inverting KG operator
void multiply_by_x_space_Wilson_gluon_propagator_by_inv(spin1prop *prop_out,spin1prop *prop_in,gluon_info gl)
{
  //source and temp prop
  spin1field *source=nissa_malloc("source",loc_vol+bord_vol,spin1field);
  spin1field *tprop=nissa_malloc("tprop",loc_vol+bord_vol,spin1field);
  
  //loop over the source index
  for(int mu_so=0;mu_so<4;mu_so++)
    {
      //copy the in
      NISSA_LOC_VOL_LOOP(ivol)
        for(int mu_si=0;mu_si<4;mu_si++)
          memcpy(source[ivol][mu_si],prop_in[ivol][mu_si][mu_so],sizeof(complex));
      set_borders_invalid(source);
      
      //invert
      inv_Wilson_gluon_Klein_Gordon_operator(tprop,NULL,gl,1000000,5,1.e-26,source);
      
      //copy into the out
      NISSA_LOC_VOL_LOOP(ivol)
        for(int mu_si=0;mu_si<4;mu_si++)
          memcpy(prop_out[ivol][mu_si][mu_so],tprop[ivol][mu_si],sizeof(complex));
    }
  
  set_borders_invalid(prop_out);
  
  nissa_free(source);
  nissa_free(tprop);
}

//compute the Wilson action gluon propagator in the x space by inverting KG operator
void compute_x_space_Wilson_gluon_propagator_by_inv(spin1prop *prop,gluon_info gl)
{
  spin1prop *source=nissa_malloc("source",loc_vol+bord_vol,spin1prop);
  
  //prepare the source
  memset(source,0,sizeof(spin1prop)*loc_vol);
  if(rank==0) spinspin_put_to_id(source[0]);
  
  //multiply
  multiply_by_x_space_Wilson_gluon_propagator_by_inv(prop,source,gl);
  
  nissa_free(source);
}
