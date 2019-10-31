#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

#include "fft.hpp"

namespace nissa
{
  //interpret free index as source or sink
  THREADABLE_FUNCTION_6ARG(pass_spinspin_from_mom_to_x_space, spinspin*,out, spinspin*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    //compute the main part of the fft
    //+1 if sink, -1 if source
    int sign[2]={-1,+1};
    int s=sign[source_or_sink]*include_phases;
    fft4d((complex*)out,(complex*)in,dirs,sizeof(spinspin)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	safe_spinspin_prod_complex(out[ivol],out[ivol],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //interprets free index as source or sink (see above)
  THREADABLE_FUNCTION_6ARG(pass_spinspin_from_x_to_mom_space, spinspin*,out, spinspin*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    //-1 if sink, +1 if source
    int sign[2]={+1,-1};
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase and put 1/vol
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	safe_spinspin_prod_complex(out[ivol],in[ivol],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spinspin)/sizeof(complex),sign[source_or_sink],1);
  }
  THREADABLE_FUNCTION_END
  
  //see above
  THREADABLE_FUNCTION_6ARG(pass_spin1prop_from_mom_to_x_space, spin1prop*,out, spin1prop*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    //+1 if sink, -1 if source
    int sign[2]={-1,+1};
    int s=sign[source_or_sink]*include_phases;
    
    //multiply by exp(i sign*(p_mu-p_nu)/2)
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
      {
	complex ph[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    double pmu=dirs[mu]*sign[source_or_sink]*M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu]*include_phases)/glb_size[mu];
	    double pmuh=pmu*0.5;
	    ph[mu][RE]=cos(pmuh);
	    ph[mu][IM]=sin(pmuh);
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  for(int nu=0;nu<NDIM;nu++)
	    {
	      safe_complex_prod      (out[imom][mu][nu],in[imom][mu][nu],ph[mu]);
	      safe_complex_conj2_prod(out[imom][mu][nu],out[imom][mu][nu],ph[nu]);
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1prop)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=dirs[mu]*steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int mu=0;mu<NDIM;mu++)
	  for(int nu=0;nu<NDIM;nu++)
	    safe_complex_prod(out[ivol][mu][nu],out[ivol][mu][nu],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //see previous note
  THREADABLE_FUNCTION_6ARG(pass_spin1prop_from_x_to_mom_space, spin1prop*,out, spin1prop*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    //-1 if sink, +1 if source
    int sign[2]={+1,-1};
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=s*dirs[mu]*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int mu=0;mu<NDIM;mu++)
	  for(int nu=0;nu<NDIM;nu++)
	    safe_complex_prod(out[ivol][mu][nu],in[ivol][mu][nu],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1prop)/sizeof(complex),sign[source_or_sink],1);
    
    //multiply by exp(i -(p_mu-p_nu)/2) and put 1/vol
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
      {
	complex ph[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    double pmu=sign[source_or_sink]*dirs[mu]*M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu]*include_phases)/glb_size[mu];
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
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_6ARG(pass_spin1field_from_mom_to_x_space, spin1field*,out, spin1field*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    //+1 if sink, -1 if source
    int sign[2]={-1,+1};
    int s=sign[source_or_sink]*include_phases;
    
    //multiply by exp(i p_mu/2)
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
      {
	complex ph[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    double pmu=dirs[mu]*sign[source_or_sink]*M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu]*include_phases)/glb_size[mu];
	    double pmuh=pmu*0.5;
	    ph[mu][RE]=cos(pmuh);
	    ph[mu][IM]=sin(pmuh);
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[imom][mu],in[imom][mu],ph[mu]);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1field)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[ivol][mu],out[ivol][mu],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_6ARG(pass_spin1field_from_x_to_mom_space, spin1field*,out, spin1field*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    //-1 if sink, +1 if source
    int sign[2]={+1,-1};
    int s=sign[source_or_sink]*include_phases;
	
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[ivol][mu],in[ivol][mu],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1field)/sizeof(complex),sign[source_or_sink],1);
    
    //multiply by exp(-i p_mu/2)
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
      {
	complex ph[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    double pmu=dirs[mu]*sign[source_or_sink]*M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu]*include_phases)/glb_size[mu];
	    double pmuh=pmu*0.5;
	    ph[mu][RE]=cos(pmuh);
	    ph[mu][IM]=sin(pmuh);
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[imom][mu],out[imom][mu],ph[mu]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_6ARG(pass_spin_from_mom_to_x_space, spin*,out, spin*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    int sign[2]={-1,+1};
    int s=sign[source_or_sink]*include_phases;
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,dirs,sizeof(spin)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  safe_complex_prod(out[ivol][id],out[ivol][id],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_6ARG(pass_spin_from_x_to_mom_space, spin*,out, spin*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    int sign[2]={+1,-1};
    int s=sign[source_or_sink]*include_phases;
       
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  safe_complex_prod(out[ivol][id],in[ivol][id],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin)/sizeof(complex),sign[source_or_sink],1);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_6ARG(pass_spincolor_from_mom_to_x_space, spincolor*,out, spincolor*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    int sign[2]={-1,+1};
    int s=sign[source_or_sink]*include_phases;
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,dirs,sizeof(spincolor)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    safe_complex_prod(out[ivol][id][ic],out[ivol][id][ic],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_6ARG(pass_spincolor_from_x_to_mom_space, spincolor*,out, spincolor*,in, bool*,dirs, double*,bc, int,source_or_sink, bool,include_phases)
  {
    GET_THREAD_ID();
    
    int sign[2]={+1,-1};
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    momentum_t steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glb_size[mu];
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    safe_complex_prod(out[ivol][id][ic],in[ivol][id][ic],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spincolor)/sizeof(complex),sign[source_or_sink],1);
  }
  THREADABLE_FUNCTION_END
}
