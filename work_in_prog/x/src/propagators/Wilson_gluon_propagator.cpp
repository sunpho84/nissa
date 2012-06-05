#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.h"
#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/complex.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/operations/fft.h"
#include "../../../../src/base/vectors.h"

#include "../types/types.h"
#include "../routines/fourier.h"
#include "../inverters/cg_Wilson_gluon_operator.h"

//compute the Wilson action gluon propagator in the momentum space according to P.Weisz
void compute_mom_space_Wilson_gluon_propagator(spin1prop *prop,gluon_info gl)
{
  //check absence of zero modes
  int zmpres=1;
  for(int mu=0;mu<4;mu++) zmpres&=(gl.bc[mu]==0);
  //if(zmpres) crash("zero mode present, prop not defined");
  
  //reset the propagator
  memset(prop,0,loc_vol*sizeof(spin1prop));
  
  nissa_loc_vol_loop(imom)
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

//compute the Wilson action gluon propagator in the x space by inverting KG operator
void compute_x_space_Wilson_gluon_propagator_by_inv(spin1prop *prop,gluon_info gl)
{
  //source and temp prop
  spin1field *source=nissa_malloc("source",loc_vol+bord_vol,spin1field);
  spin1field *tprop=nissa_malloc("tprop",loc_vol+bord_vol,spin1field);
  
  //loop over the source index
  for(int mu_so=0;mu_so<4;mu_so++)
    {
      //prepare the source
      memset(source,0,sizeof(spin1field)*loc_vol);
      if(rank==0) source[0][mu_so][0]=1;
      set_borders_invalid(source);

      //invert and copy into the spinspin
      inv_Wilson_gluon_Klein_Gordon_operator(tprop,NULL,gl,1000000,5,1.e-26,source);
      nissa_loc_vol_loop(ivol)
        for(int mu_si=0;mu_si<4;mu_si++)
          memcpy(prop[ivol][mu_si][mu_so],tprop[ivol][mu_si],sizeof(complex));
    }
  
  set_borders_invalid(prop);
  
  nissa_free(source);
  nissa_free(tprop);

}
