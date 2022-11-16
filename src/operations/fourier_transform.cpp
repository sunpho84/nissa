#ifdef HAVE_CONFIG_H
# include "config.hpp"
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
  void pass_spinspin_from_mom_to_x_space(spinspin* out,spinspin* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    //compute the main part of the fft
    //+1 if sink, -1 if source
    int sign[2]={-1,+1};
    int s=sign[source_or_sink]*include_phases;
    fft4d((complex*)out,(complex*)in,dirs,sizeof(spinspin)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	safe_spinspin_prod_complex(out[ivol.nastyConvert()],out[ivol.nastyConvert()],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //interprets free index as source or sink (see above)
  void pass_spinspin_from_x_to_mom_space(spinspin* out,spinspin* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    //-1 if sink, +1 if source
    int sign[2]={+1,-1};
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase and put 1/vol
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	safe_spinspin_prod_complex(out[ivol.nastyConvert()],in[ivol.nastyConvert()],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spinspin)/sizeof(complex),sign[source_or_sink],1);
  }
  
  //see above
  void pass_spin1prop_from_mom_to_x_space(spin1prop* out,spin1prop* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    //+1 if sink, -1 if source
    int _sign[2]={-1,+1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //multiply by exp(i sign*(p_mu-p_nu)/2)
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	complex ph[NDIM];
	FOR_ALL_DIRS(mu)
	  {
	    double pmu=dirs(mu)*sign[source_or_sink]*M_PI*(2*glbCoordOfLoclx(imom,mu)()+bc(mu)*include_phases)/glbSize(mu)();
	    double pmuh=pmu*0.5;
	    ph[mu.nastyConvert()][RE]=cos(pmuh);
	    ph[mu.nastyConvert()][IM]=sin(pmuh);
	  }
	
	FOR_ALL_DIRS(mu)
	  FOR_ALL_DIRS(nu)
	    {
	      safe_complex_prod      (out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],in[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],ph[mu.nastyConvert()]);
	      safe_complex_conj2_prod(out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],ph[nu.nastyConvert()]);
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1prop)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=dirs(mu)*steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	FOR_ALL_DIRS(mu)
	  FOR_ALL_DIRS(nu)
	    safe_complex_prod(out[ivol.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],out[ivol.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //see previous note
  void pass_spin1prop_from_x_to_mom_space(spin1prop* out,spin1prop* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    //-1 if sink, +1 if source
    int _sign[2]={+1,-1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=s*dirs(mu)*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	FOR_ALL_DIRS(mu)
	  FOR_ALL_DIRS(nu)
	    safe_complex_prod(out[ivol.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],in[ivol.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1prop)/sizeof(complex),sign[source_or_sink],1);
    
    //multiply by exp(i -(p_mu-p_nu)/2) and put 1/vol
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	complex ph[NDIM];
	FOR_ALL_DIRS(mu)
	  {
	    const double pmu=sign[source_or_sink]*dirs(mu)*M_PI*(2*glbCoordOfLoclx(imom,mu)()+bc(mu)*include_phases)/glbSize(mu)();
	    const double pmuh=pmu*0.5;
	    ph[mu.nastyConvert()][RE]=cos(pmuh);
	    ph[mu.nastyConvert()][IM]=sin(pmuh);
	  }
	
	FOR_ALL_DIRS(mu)
	  FOR_ALL_DIRS(nu)
	    {
	      safe_complex_prod      (out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],ph[mu.nastyConvert()]);
	      safe_complex_conj2_prod(out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],out[imom.nastyConvert()][mu.nastyConvert()][nu.nastyConvert()],ph[nu.nastyConvert()]);
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void pass_spin1field_from_mom_to_x_space(spin1field* out,spin1field* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    //+1 if sink, -1 if source
    int _sign[2]={-1,+1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //multiply by exp(i p_mu/2)
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	complex ph[NDIM];
	FOR_ALL_DIRS(mu)
	  {
	    const double pmu=dirs(mu)*sign[source_or_sink]*M_PI*(2*glbCoordOfLoclx(imom,mu)()+bc(mu)*include_phases)/glbSize(mu)();
	    const double pmuh=pmu*0.5;
	    ph[mu.nastyConvert()][RE]=cos(pmuh);
	    ph[mu.nastyConvert()][IM]=sin(pmuh);
	  }
	
	FOR_ALL_DIRS(mu)
	  safe_complex_prod(out[imom.nastyConvert()][mu.nastyConvert()],in[imom.nastyConvert()][mu.nastyConvert()],ph[mu.nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1field)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	FOR_ALL_DIRS(mu)
	  safe_complex_prod(out[ivol.nastyConvert()][mu.nastyConvert()],out[ivol.nastyConvert()][mu.nastyConvert()],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void pass_spin1field_from_x_to_mom_space(spin1field* out,spin1field* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    //-1 if sink, +1 if source
    int _sign[2]={+1,-1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	FOR_ALL_DIRS(mu)
	  safe_complex_prod(out[ivol.nastyConvert()][mu.nastyConvert()],in[ivol.nastyConvert()][mu.nastyConvert()],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1field)/sizeof(complex),sign[source_or_sink],1);
    
    //multiply by exp(-i p_mu/2)
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	complex ph[NDIM];
	FOR_ALL_DIRS(mu)
	  {
	    const double pmu=dirs(mu)*sign[source_or_sink]*M_PI*(2*glbCoordOfLoclx(imom,mu)()+bc(mu)*include_phases)/glbSize(mu)();
	    double pmuh=pmu*0.5;
	    ph[mu.nastyConvert()][RE]=cos(pmuh);
	    ph[mu.nastyConvert()][IM]=sin(pmuh);
	  }
	
	FOR_ALL_DIRS(mu)
	  safe_complex_prod(out[imom.nastyConvert()][mu.nastyConvert()],out[imom.nastyConvert()][mu.nastyConvert()],ph[mu.nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void pass_spin_from_mom_to_x_space(spin* out,spin* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    int _sign[2]={-1,+1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,dirs,sizeof(spin)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  safe_complex_prod(out[ivol.nastyConvert()][id],out[ivol.nastyConvert()][id],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void pass_spin_from_x_to_mom_space(spin* out,spin* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    int _sign[2]={+1,-1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
       
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  safe_complex_prod(out[ivol.nastyConvert()][id],in[ivol.nastyConvert()][id],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spin)/sizeof(complex),sign[source_or_sink],1);
  }
  
  void pass_spincolor_from_mom_to_x_space(spincolor* out,spincolor* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    int _sign[2]={-1,+1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)in,dirs,sizeof(spincolor)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    safe_complex_prod(out[ivol.nastyConvert()][id][ic],out[ivol.nastyConvert()][id][ic],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  void pass_spincolor_from_x_to_mom_space(spincolor* out,spincolor* in,const Coords<bool>& dirs,const Momentum& bc,int source_or_sink,bool include_phases)
  {
    
    int _sign[2]={+1,-1},*sign=_sign;
    int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    Momentum steps;
    FOR_ALL_DIRS(mu)
      steps(mu)=dirs(mu)*s*bc(mu)*M_PI/glbSize(mu)();
    
    //add the fractional phase
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute phase exponent
	double arg=0;
	FOR_ALL_DIRS(mu)
	  arg+=steps(mu)*glbCoordOfLoclx(ivol,mu)();
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    safe_complex_prod(out[ivol.nastyConvert()][id][ic],in[ivol.nastyConvert()][id][ic],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //compute the main part of the fft
    fft4d((complex*)out,(complex*)out,dirs,sizeof(spincolor)/sizeof(complex),sign[source_or_sink],1);
  }
}
