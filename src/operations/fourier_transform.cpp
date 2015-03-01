#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/macros.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"

#include "fft.hpp"

namespace nissa
{
  void pass_spinspin_from_mom_to_x_space(spinspin *out,spinspin *in,momentum_t bc)
  {
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,16,+1,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      safe_spinspin_complex_prod(out[ivol],out[ivol],ph);
    }
  }
  
  void pass_spinspin_from_x_to_mom_space(spinspin *out,spinspin *in,momentum_t bc)
  {
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=-bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++)
	arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
      
      //compute the phase and put 1/vol
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      safe_spinspin_complex_prod(out[ivol],in[ivol],ph);
    }
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,16,-1,1);  
  }
  
  void pass_spin1prop_from_mom_to_x_space(spin1prop *out,spin1prop *in,momentum_t bc)
  {
    //multiply by exp(i (p_mu-p_nu)/2)
    NISSA_LOC_VOL_LOOP(imom)
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
	    safe_complex_conj2_prod(out[imom][mu][nu],out[imom][mu][nu],ph[nu]);
	  }
    }
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,16,+1,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
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
    NISSA_LOC_VOL_LOOP(ivol)
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
    fft4d((complex*)out,(complex*)out,16,-1,1);
    
    //multiply by exp(i -(p_mu-p_nu)/2) and put 1/vol
    NISSA_LOC_VOL_LOOP(imom)
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
  
  void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,momentum_t bc,bool bar=false)
  {
    int sign=+1;
    if(bar) sign*=-1;
    
    //multiply by exp(i p_mu/2)
    NISSA_LOC_VOL_LOOP(imom)
    {
      complex ph[4];
      for(int mu=0;mu<4;mu++)
	{
	  double pmu=sign*M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  double pmuh=pmu*0.5;
	  ph[mu][RE]=cos(pmuh);
	  ph[mu][IM]=sin(pmuh);
	}
      
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[imom][mu],in[imom][mu],ph[mu]);
    }
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,4,sign,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
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
  
  void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,momentum_t bc,bool bar=false)
  {
    int sign=-1;
    if(bar) sign*=-1;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
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
    fft4d((complex*)out,(complex*)out,4,sign,1);
    
    //multiply by exp(-i p_mu/2)
    NISSA_LOC_VOL_LOOP(imom)
    {
      complex ph[4];
      for(int mu=0;mu<4;mu++)
	{
	  double pmu=sign*M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  double pmuh=pmu*0.5;
	  ph[mu][RE]=cos(pmuh);
	  ph[mu][IM]=sin(pmuh);
	}
      
      for(int mu=0;mu<4;mu++)
	safe_complex_prod(out[imom][mu],out[imom][mu],ph[mu]);
    }
  }
  
  void pass_spin_from_mom_to_x_space(spin *out,spin *in,momentum_t bc,bool bar=false)
  {
    int sign=+1;
    if(bar) sign*=-1;
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,4,+sign,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
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
  
  void pass_spin_from_x_to_mom_space(spin *out,spin *in,momentum_t bc,bool bar=false)
  {
    int sign=-1;
    if(bar) sign*=-1;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_LOC_VOL_LOOP(ivol)
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
    fft4d((complex*)out,(complex*)out,4,+sign,1);
  }
}
