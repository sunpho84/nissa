#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.h"
#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/complex.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/operations/fft.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/base/communicate.h"

#include "../types/types.h"
#include "../routines/fourier.h"

void apply_Wilson_gluon_x_Klein_Gordon_operator(spin1field *out,spin1field *in,gluon_info gl)
{
  communicate_lx_spin_borders((spin*)in);
  
  double rep_alpha=1/gl.alpha;
  
  //this is the field where we will store the backward derivative of the "in" field
  //the derivative is taken with respect to the first index
  spin1field *bder[4];
  for(int mu=0;mu<4;mu++) bder[mu]=nissa_malloc("bder",loc_vol+bord_vol,spin1field);
  
  //compute the border phases
  complex phases[4];
  for(int mu=0;mu<4;mu++)
    {
      double arg=M_PI*gl.bc[mu];
      phases[mu][RE]=cos(arg);
      phases[mu][IM]=sin(arg);
    }

  //backward derivative
  nissa_loc_vol_loop(ivol)
    //direction of the derivative
    for(int mu=0;mu<4;mu++)
      {
	//neighbour
	int idw=loclx_neighdw[ivol][mu];
	
	//check if we are on lower border of direction mu
	int lbord=(glb_coord_of_loclx[ivol][mu]==0);
	
	//compute factor
	complex fact={1,0};
	if(lbord) complex_copy(fact,phases[mu]);
	
	//loop on the four indices
	for(int nu=0;nu<4;nu++)
	  {
	    complex_copy               (bder[mu][ivol][nu],in[ivol][nu]);
	    complex_subt_the_conj2_prod(bder[mu][ivol][nu], in[idw][nu],fact);
	  }
      }

  for(int mu=0;mu<4;mu++) communicate_lx_spin_borders(bder[mu]);
  
  nissa_loc_vol_loop(ivol)
    //index of the field
    for(int mu=0;mu<4;mu++)
      {
	//reset the output
	memset(out[ivol][mu],0,sizeof(complex));
	
	for(int nu=0;nu<4;nu++)
	  {
	    //part1: forward derivative on the back derivative direction
	    //derivative index: nu
	    int iup=loclx_neighup[ivol][nu];

	    //check if we are on upper border of direction nu
	    int ubord=(glb_coord_of_loclx[ivol][nu]==glb_size[nu]-1);
	    complex fact1={1,0};
	    if(ubord) memcpy(fact1,phases[nu],sizeof(complex));

	    complex_summ_the_prod(out[ivol][mu],bder[nu] [iup][mu],fact1);
	    complex_subtassign   (out[ivol][mu],bder[nu][ivol][mu]);
	    
	    //part2: forward derivative on field index direction
	    //derivative index: mu
	    iup=loclx_neighup[ivol][mu];
	    
	    //check if we are on upper border of direction mu
	    ubord=(glb_coord_of_loclx[iup][mu]==0);
	    double coef2=rep_alpha-1;
	    complex fact2={coef2,0};
	    if(ubord) safe_complex_prod(fact2,fact2,phases[mu]);
	    complex_summ_the_prod       (out[ivol][mu],bder[nu] [iup][nu],fact2);
	    complex_subt_the_prod_double(out[ivol][mu],bder[nu][ivol][nu],coef2);
	  }
      }
  
  //put minus sign
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      complex_prodassign_double(out[ivol][mu],-1);
  
  for(int mu=0;mu<4;mu++) nissa_free(bder[mu]);
}

void apply_Wilson_gluon_mom_Klein_Gordon_operator(spin1field *out,spin1field *in,gluon_info gl)
{
  double rep_alpha=1/gl.alpha;
  int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  
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
}
