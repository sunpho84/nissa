#include <string.h>
#include <math.h>

#include "../../src/base/global_variables.h"
#include "../../src/new_types/new_types_definitions.h"
#include "../../src/new_types/complex.h"
#include "../../src/base/debug.h"
#include "../../src/operations/fft.h"
#include "types.h"

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
  
  //multiply by exp(i\sum_mu p_mu/2)
  nissa_loc_vol_loop(imom)
    {
      double s=0;
      for(int mu=0;mu<4;mu++) s+=glb_coord_of_loclx[imom][mu];
      s/=2;
      complex ph={cos(s),sin(s)};
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(prop[imom][id1][id2],prop[imom][id1][id2],ph);
    }
  
  //compute the main part of the fft
  fft4d((complex*)prop,(complex*)prop,16,+1,1);
  
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++) steps[mu]=gl.bc[mu]*M_PI/glb_size[mu];

  //add the fractional phase
  nissa_loc_vol_loop(imom)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++) arg+=steps[mu]*(glb_coord_of_loclx[imom][mu]+0.5);
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(prop[imom][id1][id2],prop[imom][id1][id2],ph);
    }
}
