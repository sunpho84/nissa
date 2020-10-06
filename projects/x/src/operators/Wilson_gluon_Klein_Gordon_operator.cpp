#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.hpp"
#include "../../../../src/new_types/new_types_definitions.hpp"
#include "../../../../src/new_types/complex.hpp"
#include "../../../../src/base/debug.hpp"
#include "../../../../src/operations/fft.hpp"
#include "../../../../src/base/vectors.hpp"
#include "../../../../src/communicate/communicate.hpp"

#include "../types/types.hpp"
#include "../routines/derivatives.hpp"

void apply_Wilson_gluon_x_Klein_Gordon_operator(spin1field *out,spin1field *in,gluon_info gl)
{
  if(out==in) crash("in==out");
  
  double rep_alpha=1/gl.alpha;
  
  //this is the field where we will store the backward derivative of the "in" field
  //the derivative is taken with respect to the first index
  spin1field *bder=nissa_malloc("bder",loc_vol+bord_vol,spin1field),*fbder[4][4];
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      fbder[mu][nu]=nissa_malloc("fbder",loc_vol,spin1field);
  
  //reset the output
  memset(out,0,sizeof(spin1field)*loc_vol);
  
  //second order derivative
  for(int nu=0;nu<4;nu++)
    {
      spin1field_bw_derivative(bder,in,gl.bc,nu);
      for(int mu=0;mu<4;mu++) spin1field_fw_derivative(fbder[mu][nu],bder,gl.bc,mu);
    }

  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	{
	  //part1: forward derivative on the back derivative direction
	  //forward derivative index: nu, backward: nu, field: mu
	  complex_summassign(out[ivol][mu],fbder[nu][nu] [ivol] [mu]);
	  //part2: forward derivative on field index direction
	  //forward derivative index: mu, backward: nu, field: nu
	  complex_summ_the_prod_double(out[ivol][mu],fbder[mu][nu] [ivol] [nu],rep_alpha-1);
	}
  
  //put minus sign
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      complex_prodassign_double(out[ivol][mu],-1);
  
  set_borders_invalid(out);
  
  for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) nissa_free(fbder[mu][nu]);
  nissa_free(bder);
}

void apply_Wilson_gluon_mom_Klein_Gordon_operator(spin1field *out,spin1field *in,gluon_info gl)
{
  double rep_alpha=1/gl.alpha;
  int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  
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
      
      spin1field temp;
      for(int mu=0;mu<4;mu++)
	for(int ri=0;ri<2;ri++)
	  {
	    temp[mu][ri]=0;
	    for(int nu=0;nu<4;nu++)
	      temp[mu][ri]+=(kron_delta[mu][nu]*kt2+(rep_alpha-1)*kt[mu]*kt[nu])*in[imom][nu][ri];
	  }
      memcpy(out[imom],temp,sizeof(spin1field));
    }
  
  set_borders_invalid(out);
}
