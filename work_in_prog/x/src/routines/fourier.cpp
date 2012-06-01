#include <math.h>

#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/complex.h"
#include "../../../../src/operations/fft.h"
#include "../../../../src/base/global_variables.h"

void pass_spinspin_from_mom_to_x_space(spinspin *out,spinspin *in,momentum_t bc)
{
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)in,16,+1,1);
  
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(out[ivol][id1][id2],out[ivol][id1][id2],ph);
    }
}

void pass_spinspin_from_x_to_mom_space(spinspin *out,spinspin *in,momentum_t bc)
{
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=-bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(out[ivol][id1][id2],in[ivol][id1][id2],ph);
    }
  
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)out,16,-1,0);  
}

void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,momentum_t bc)
{
  //multiply by exp(i (p_mu-p_nu)/2)
  nissa_loc_vol_loop(imom)
    {
      complex ph[4];
      for(int mu=0;mu<4;mu++)
	{
	  double pmu=M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  double pmuh=pmu*0.5;
	  ph[mu][RE]=cos(pmuh);
	  ph[mu][IM]=sin(pmuh);
	}
      
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    safe_complex_prod      (out[imom][mu][nu],in[imom][mu][nu],ph[mu]);
	    safe_complex_conj2_prod(out[imom][mu][nu],in[imom][mu][nu],ph[nu]);
	  }
    }
  
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)out,16,+1,1);
  
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  safe_complex_prod(out[ivol][mu][nu],out[ivol][mu][nu],ph);
    }
}

void pass_spin1prop_from_x_to_mom_space(spin1prop *out,spin1prop *in,momentum_t bc)
{
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=-bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  safe_complex_prod(out[ivol][mu][nu],in[ivol][mu][nu],ph);
    }
  
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)out,16,-1,0);
  
  //multiply by exp(i -(p_mu-p_nu)/2)
  nissa_loc_vol_loop(imom)
    {
      complex ph[4];
      for(int mu=0;mu<4;mu++)
	{
	  double pmu=-M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  double pmuh=pmu*0.5;
	  ph[mu][RE]=cos(pmuh);
	  ph[mu][IM]=sin(pmuh);
	}
      
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    safe_complex_prod      (out[imom][mu][nu],out[imom][mu][nu],ph[mu]);
	    safe_complex_conj2_prod(out[imom][mu][nu],out[imom][mu][nu],ph[nu]);
	  }
    }
}

void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,momentum_t bc)
{
  //multiply by exp(i (p_mu-p_nu)/2)
  nissa_loc_vol_loop(imom)
    {
      complex ph[4];
      for(int mu=0;mu<4;mu++)
	{
	  double pmu=M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  double pmuh=pmu*0.5;
	  ph[mu][RE]=cos(pmuh);
	  ph[mu][IM]=sin(pmuh);
	}
      
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[imom][mu],in[imom][mu],ph[mu]);
    }
  
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)out,4,+1,1);
  
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[ivol][mu],out[ivol][mu],ph);
    }
}

void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,momentum_t bc)
{
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=-bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[ivol][mu],in[ivol][mu],ph);
    }
  
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)out,4,-1,0);
  
  //multiply by exp(-i p_mu/2)
  nissa_loc_vol_loop(imom)
    {
      complex ph[4];
      for(int mu=0;mu<4;mu++)
	{
	  double pmu=-M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  double pmuh=pmu*0.5;
	  ph[mu][RE]=cos(pmuh);
	  ph[mu][IM]=sin(pmuh);
	}
      
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[imom][mu],out[imom][mu],ph[mu]);
    }
}

void pass_spin_from_mom_to_x_space(spin *out,spin *in,momentum_t bc)
{
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)in,4,+1,1);
  
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[ivol][mu],out[ivol][mu],ph);
    }
}

void pass_spin_from_x_to_mom_space(spin *out,spin *in,momentum_t bc)
{
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++)
    steps[mu]=-bc[mu]*M_PI/glb_size[mu];
  
  //add the fractional phase
  nissa_loc_vol_loop(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[ivol][mu],in[ivol][mu],ph);
    }
  
  //compute the main part of the fft
  fft4d((complex*)out,(complex*)out,4,-1,0);
}
