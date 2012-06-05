#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.h"
#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/complex.h"
#include "../../../../src/new_types/spin.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/base/routines.h"
#include "../../../../src/operations/fft.h"

#include "../types/types.h"
#include "../routines/fourier.h"

//compute the tree level Symanzik gluon propagator in the momentum space according to P.Weisz
void mom_space_tlSym_gluon_propagator_of_imom(spin1prop prop,gluon_info gl,int imom)
{
  double c1=gl.c1;
  int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  
  //momentum
  momentum_t k,kt;
  double kt2=0,kt4=0,kt6=0;
  double ktpo2[4][4],ktso2[4][4];
  for(int mu=0;mu<4;mu++)
    {
      k[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+gl.bc[mu])/glb_size[mu];
      kt[mu]=2*sin(k[mu]/2);
      kt2+=kt[mu]*kt[mu];
      kt4+=kt[mu]*kt[mu]*kt[mu]*kt[mu];
      kt6+=kt[mu]*kt[mu]*kt[mu]*kt[mu]*kt[mu]*kt[mu];
      for(int nu=0;nu<4;nu++)
	{
	  ktpo2[mu][nu]=1;
	  ktso2[mu][nu]=0;
	  for(int rho=0;rho<4;rho++)
	    if(mu!=rho && nu!=rho)
	      {
		ktpo2[mu][nu]*=kt[rho]*kt[rho];
		ktso2[mu][nu]+=kt[rho]*kt[rho];
	      }
	}
    }
  double kt22=kt2*kt2;
  double kt23=kt2*kt2*kt2;
  double kt42=kt4*kt4;
  
  if(kt2!=0)
    {
      //Deltakt
      double Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c1*c1*(kt23+2*kt6-kt2*kt4));
      for(int rho=0;rho<4;rho++) Deltakt-=4*c1*c1*c1*kt[rho]*kt[rho]*kt[rho]*kt[rho]*ktpo2[rho][rho];
      
      //A
      double A[4][4];
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  A[mu][nu]=(1-kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2*kt4+kt2*ktso2[mu][nu])+c1*c1*(kt42+kt2*kt4*ktpo2[mu][nu]+kt22*ktpo2[mu][nu]));
      
      //Prop
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    prop[mu][nu][RE]=gl.alpha*kt[mu]*kt[nu];
	    for(int si=0;si<4;si++)
	      prop[mu][nu][RE]+=(kt[si]*kron_delta[mu][nu]-kt[nu]*kron_delta[mu][si])*kt[si]*A[si][nu];
	    
	    prop[mu][nu][RE]/=kt2*kt2*glb_vol;
	    prop[mu][nu][IM]=0;
	  }
    }
  else
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	prop[mu][nu][RE]=prop[mu][nu][IM]=gl.zmp/glb_vol;
}

void compute_mom_space_tlSym_gluon_propagator(spin1prop *prop,gluon_info gl)
{
  nissa_loc_vol_loop(imom)
    mom_space_tlSym_gluon_propagator_of_imom(prop[imom],gl,imom);
  set_borders_invalid(prop);
}

void multiply_mom_space_tlSym_gluon_propagator(spin1field *out,spin1field *in,gluon_info gl)
{
  nissa_loc_vol_loop(imom)
    {
      spin1prop prop;
      mom_space_tlSym_gluon_propagator_of_imom(prop,gl,imom);
      safe_spinspin_spin_prod(out[imom],prop,in[imom]);
    }
  set_borders_invalid(out);
}

void multiply_x_space_tlSym_gluon_propagator_by_fft(spin1prop *out,spin1prop *in,gluon_info gl)
{
  pass_spin1prop_from_x_to_mom_space(out,in,gl.bc);
  nissa_loc_vol_loop(imom)
    {
      spin1prop prop;
      mom_space_tlSym_gluon_propagator_of_imom(prop,gl,imom);
      safe_spinspin_spinspin_prod(out[imom],prop,out[imom]);
    }
  pass_spin1prop_from_mom_to_x_space(out,in,gl.bc);
  set_borders_invalid(out);
}

//compute the tree level Symanzik gluon propagator in the x space by taking the fft of that in momentum space
void compute_x_space_tlSym_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl)
{
  compute_mom_space_tlSym_gluon_propagator(prop,gl);
  pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc);
}
