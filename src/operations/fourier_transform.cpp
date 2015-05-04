#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/macros.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"

#include "fft.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  THREADABLE_FUNCTION_3ARG(pass_spinspin_from_mom_to_x_space, spinspin*,out, spinspin*,in, double*,bc)
  {
    GET_THREAD_ID();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,16,+1,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_3ARG(pass_spinspin_from_x_to_mom_space, spinspin*,out, spinspin*,in, double*,bc)
  {
    GET_THREAD_ID();
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=-bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,16,-1,1);  
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_3ARG(pass_spin1prop_from_mom_to_x_space, spin1prop*,out, spin1prop*,in, double*,bc)
  {
    GET_THREAD_ID();
    
    //multiply by exp(i (p_mu-p_nu)/2)
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
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
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,16,+1,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_3ARG(pass_spin1prop_from_x_to_mom_space, spin1prop*,out, spin1prop*,in, double*,bc)
  {
    GET_THREAD_ID();
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=-bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,16,-1,1);
    
    //multiply by exp(i -(p_mu-p_nu)/2) and put 1/vol
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
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
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  void pass_spin1field_from_mom_to_x_space(spin1field *out,spin1field *in,double *bc,bool bar=false)
  {
    GET_THREAD_ID();
    
    int sign=+1;
    if(bar) sign*=-1;
    
    //multiply by exp(i p_mu/2)
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
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
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,4,sign,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    set_borders_invalid(out);
  }

  void pass_spin1field_from_x_to_mom_space(spin1field *out,spin1field *in,double *bc,bool bar=false)
  {
    GET_THREAD_ID();
    
    int sign=-1;
    if(bar) sign*=-1;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,4,sign,1);
    
    //multiply by exp(-i p_mu/2)
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
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
    set_borders_invalid(out);
  }
  
  void pass_spin_from_mom_to_x_space(spin *out,spin *in,double *bc,bool bar=false)
  {
    GET_THREAD_ID();
    
    int sign=+1;
    if(bar) sign*=-1;
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,4,+sign,0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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

    set_borders_invalid(out);
  }
  
  void pass_spin_from_x_to_mom_space(spin *out,spin *in,double *bc,bool bar=false)
  {
    GET_THREAD_ID();
    
    int sign=-1;
    if(bar) sign*=-1;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<4;mu++)
      steps[mu]=sign*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
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
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,4,+sign,1);
  }
}
