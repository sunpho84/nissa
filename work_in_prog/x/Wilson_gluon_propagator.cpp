#include <string.h>
#include <math.h>

#include "../../src/base/global_variables.h"
#include "../../src/new_types/new_types_definitions.h"
#include "../../src/new_types/complex.h"
#include "../../src/base/debug.h"
#include "../../src/operations/fft.h"

#include "types.h"
#include "fourier.h"

//compute the Wilson action gluon propagator in the momentum space according to P.Weisz
void compute_mom_space_Wilson_gluon_propagator(spin1prop *prop,gluon_info gl)
{
  //check absence of zero modes
  int zmpres=1;
  for(int mu=0;mu<4;mu++) zmpres&=(gl.bc[mu]==0);
  if(zmpres) crash("zero mode present, prop not defined");
  
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
	    if(mu==nu) prop[imom][mu][nu][RE]=1;
	    prop[imom][mu][nu][RE]-=(1-gl.alpha)*kt[mu]*kt[nu]/kt2;
	    prop[imom][mu][nu][RE]/=kt2;
	  }
    }
}

//compute the Wilson action gluon propagator in the x space by taking the fft of that in momentum space
void compute_x_space_Wilson_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl)
{
  compute_mom_space_Wilson_gluon_propagator(prop,gl);
  pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc);
  pass_spin1prop_from_x_to_mom_space(prop,prop,gl.bc);
  pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc);
}
