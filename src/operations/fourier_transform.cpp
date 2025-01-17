#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/field.hpp"
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
  void pass_spinspin_from_mom_to_x_space(LxField<spinspin>& out,
					 const LxField<spinspin>& in,
					 const WhichDirs& dirs,
					 const Momentum& bc,
					 const int& source_or_sink,
					 const bool& include_phases)
  {
    CRASH("reimplement");
    //compute the main part of the fft
    //+1 if sink, -1 if source
    // constexpr std::array<int,2> sign={-1,+1};
    // const int s=sign[source_or_sink]*include_phases;
    
    //fft4d((complex*)out,(complex*)in,dirs,sizeof(spinspin)/sizeof(complex),sign[source_or_sink],0);
    
    // //compute steps
    // Momentum steps;
    // for(int mu=0;mu<NDIM;mu++)
    //   steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    // //add the fractional phase
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	//compute phase exponent
    // 	double arg=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
    // 	//compute the phase
    // 	complex ph={cos(arg),sin(arg)};
	
    // 	//adapt the phase
    // 	safe_spinspin_prod_complex(out[ivol],out[ivol],ph);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
  }
  
  //interprets free index as source or sink (see above)
  void pass_spinspin_from_x_to_mom_space(LxField<spinspin>& out,
					 const LxField<spinspin>& in,
					 const WhichDirs& dirs,
					 const Momentum& bc,
					 const int& source_or_sink,
					 const bool& include_phases)
  {
    
    // //-1 if sink, +1 if source
    // constexpr std::array<int,2> sign={+1,-1};
    // const int s=sign[source_or_sink]*include_phases;
    
    // //compute steps
    // Momentum steps;
    // for(int mu=0;mu<NDIM;mu++)
    //   steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    // //add the fractional phase
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	//compute phase exponent
    // 	double arg=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
    // 	//compute the phase and put 1/vol
    // 	complex ph={cos(arg),sin(arg)};
	
    // 	//adapt the phase
    // 	safe_spinspin_prod_complex(out[ivol],in[ivol],ph);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
    
    //compute the main part of the fft
    CRASH("reimplement");
    //fft4d((complex*)out,(complex*)out,dirs,sizeof(spinspin)/sizeof(complex),sign[source_or_sink],1);
  }
  
  //see above
  void pass_spin1prop_from_mom_to_x_space(LxField<spin1prop>& out,
					  const LxField<spin1prop>& in,
					  const WhichDirs& dirs,
					  const Momentum& bc,
					  const int& source_or_sink,
					  const bool& include_phases)
  {
    // //+1 if sink, -1 if source
    // constexpr std::array<int,2> sign={-1,+1};
    // const int s=sign[source_or_sink]*include_phases;
    
    // //multiply by exp(i sign*(p_mu-p_nu)/2)
    // NISSA_PARALLEL_LOOP(imom,0,locVol)
    //   {
    // 	complex ph[NDIM];
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  {
    // 	    double pmu=dirs[mu]*sign[source_or_sink]*M_PI*(2*glbCoordOfLoclx[imom][mu]+bc[mu]*include_phases)/glbSize[mu];
    // 	    double pmuh=pmu*0.5;
    // 	    ph[mu][RE]=cos(pmuh);
    // 	    ph[mu][IM]=sin(pmuh);
    // 	  }
	
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  for(int nu=0;nu<NDIM;nu++)
    // 	    {
    // 	      safe_complex_prod      (out[imom][mu][nu],in[imom][mu][nu],ph[mu]);
    // 	      safe_complex_conj2_prod(out[imom][mu][nu],out[imom][mu][nu],ph[nu]);
    // 	    }
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
    
    //compute the main part of the fft
    CRASH("reimplement");
    //fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1prop)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    // Momentum steps;
    // for(int mu=0;mu<NDIM;mu++)
    //   steps[mu]=s*bc[mu]*M_PI/glbSize[mu];
    
    // //add the fractional phase
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	//compute phase exponent
    // 	double arg=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  arg+=dirs[mu]*steps[mu]*glbCoordOfLoclx[ivol][mu];
	
    // 	//compute the phase
    // 	complex ph={cos(arg),sin(arg)};
	
    // 	//adapt the phase
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  for(int nu=0;nu<NDIM;nu++)
    // 	    safe_complex_prod(out[ivol][mu][nu],out[ivol][mu][nu],ph);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
  }
  
  //see previous note
  void pass_spin1prop_from_x_to_mom_space(LxField<spin1prop>& out,
					  const LxField<spin1prop>& in,
					  const WhichDirs& dirs,
					  const Momentum& bc,
					  const int& source_or_sink,
					  const bool& include_phases)
  {
    
    // //-1 if sink, +1 if source
    // constexpr std::array<int,2> sign={+1,-1};
    // const int s=sign[source_or_sink]*include_phases;
    
    // //compute steps
    // Momentum steps;
    // for(int mu=0;mu<NDIM;mu++)
    //   steps[mu]=s*dirs[mu]*bc[mu]*M_PI/glbSize[mu];
    
    // //add the fractional phase
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	//compute phase exponent
    // 	double arg=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
    // 	//compute the phase
    // 	complex ph={cos(arg),sin(arg)};
	
    // 	//adapt the phase
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  for(int nu=0;nu<NDIM;nu++)
    // 	    safe_complex_prod(out[ivol][mu][nu],in[ivol][mu][nu],ph);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
    
    //compute the main part of the fft
    CRASH("reimplement");
    // fft4d((complex*)out,(complex*)out,dirs,sizeof(spin1prop)/sizeof(complex),sign[source_or_sink],1);
    
    //multiply by exp(i -(p_mu-p_nu)/2) and put 1/vol
    // NISSA_PARALLEL_LOOP(imom,0,locVol)
    //   {
    // 	complex ph[NDIM];
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  {
    // 	    double pmu=sign[source_or_sink]*dirs[mu]*M_PI*(2*glbCoordOfLoclx[imom][mu]+bc[mu]*include_phases)/glbSize[mu];
    // 	    double pmuh=pmu*0.5;
    // 	    ph[mu][RE]=cos(pmuh);
    // 	    ph[mu][IM]=sin(pmuh);
    // 	  }
	
    // 	for(int mu=0;mu<4;mu++)
    // 	  for(int nu=0;nu<4;nu++)
    // 	    {
    // 	      safe_complex_prod      (out[imom][mu][nu],out[imom][mu][nu],ph[mu]);
    // 	      safe_complex_conj2_prod(out[imom][mu][nu],out[imom][mu][nu],ph[nu]);
    // 	    }
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
  }
  
  void pass_spin1field_from_mom_to_x_space(LxField<spin1field>& out,
					   const LxField<spin1field>& in,
					   const WhichDirs& dirs,
					   const Momentum& bc,
					   const int& source_or_sink,
					   const bool& include_phases)
  {
    
    //+1 if sink, -1 if source
    constexpr std::array<int,2> sign={-1,+1};
    const int s=sign[source_or_sink]*include_phases;
    
    //multiply by exp(i p_mu/2)
    PAR(0,locVol,
	CAPTURE(dirs,include_phases,si=sign[source_or_sink],bc,
		TO_WRITE(out),
		TO_READ(in)),imom,
      {
	complex ph[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    double pmu=dirs[mu]*si*M_PI*(2*glbCoordOfLoclx[imom][mu]+bc[mu]*include_phases)/glbSize[mu];
	    double pmuh=pmu*0.5;
	    ph[mu][RE]=cos(pmuh);
	    ph[mu][IM]=sin(pmuh);
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[imom][mu],in[imom][mu],ph[mu]);
      });
    
    //compute the main part of the fft
    fft4d(out,sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    //add the fractional phase
    PAR(0,locVol,
	CAPTURE(steps,
		TO_WRITE(out)),
	ivol,
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[ivol][mu],out[ivol][mu],ph);
      });
  }
  
  void pass_spin1field_from_x_to_mom_space(LxField<spin1field>& out,
					   const LxField<spin1field>& in,
					   const WhichDirs& dirs,
					   const Momentum& bc,
					   const int& source_or_sink,
					   const bool& include_phases)
  {
    //-1 if sink, +1 if source
    constexpr std::array<int,2> sign={+1,-1};
    const int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    Momentum steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    //add the fractional phase
    PAR(0,locVol,
	CAPTURE(steps,
		TO_WRITE(out),
		TO_READ(in)),
	ivol,
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[ivol][mu],in[ivol][mu],ph);
      });
    
    //compute the main part of the fft
    fft4d(out,sign[source_or_sink],1);
    
    //multiply by exp(-i p_mu/2)
    PAR(0,locVol,
	CAPTURE(dirs,bc,include_phases,
		si=sign[source_or_sink],
		TO_WRITE(out)),
	imom,
      {
	complex ph[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    const double pmu=dirs[mu]*si*M_PI*(2*glbCoordOfLoclx[imom][mu]+bc[mu]*include_phases)/glbSize[mu];
	    const double pmuh=pmu*0.5;
	    ph[mu][RE]=cos(pmuh);
	    ph[mu][IM]=sin(pmuh);
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  safe_complex_prod(out[imom][mu],out[imom][mu],ph[mu]);
      });
  }
  
  void pass_spin_from_mom_to_x_space(LxField<spin>& out,
				     const LxField<spin>& in,
				     const WhichDirs& dirs,
				     const Momentum& bc,
				     const int& source_or_sink,
				     const bool& include_phases)
  {
    
    constexpr std::array<int,2> sign={-1,+1};
    const int s=sign[source_or_sink]*include_phases;
    
    //compute the main part of the fft
    CRASH("reimplement");
    // fft4d((complex*)out,(complex*)in,dirs,sizeof(spin)/sizeof(complex),sign[source_or_sink],0);
    
    //compute steps
    Momentum steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    //add the fractional phase
    PAR(0,locVol,
	CAPTURE(steps,
		TO_WRITE(out)),
	ivol,
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  safe_complex_prod(out[ivol][id],out[ivol][id],ph);
      });
  }
  
  void pass_spin_from_x_to_mom_space(LxField<spin>& out,
				     const LxField<spin>& in,
				     const WhichDirs& dirs,
				     const Momentum& bc,
				     const int& source_or_sink,
				     const bool& include_phases)
  {
    constexpr std::array<int,2> sign={+1,-1};
    const int s=sign[source_or_sink]*include_phases;
    
    //compute steps
    Momentum steps;
    for(int mu=0;mu<NDIM;mu++)
      steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    //add the fractional phase
    PAR(0,locVol,
	CAPTURE(steps,
		TO_WRITE(out),
		TO_READ(in)),
	ivol,
      {
	//compute phase exponent
	double arg=0;
	for(int mu=0;mu<NDIM;mu++)
	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
	//compute the phase
	complex ph={cos(arg),sin(arg)};
	
	//adapt the phase
	for(int id=0;id<NDIRAC;id++)
	  safe_complex_prod(out[ivol][id],in[ivol][id],ph);
      });
    
    //compute the main part of the fft
    CRASH("reimplement");
    // fft4d((complex*)out,(complex*)out,dirs,sizeof(spin)/sizeof(complex),sign[source_or_sink],1);
  }
  
  void pass_spincolor_from_mom_to_x_space(LxField<spincolor>& out,
					  const LxField<spincolor>& in,
					  const WhichDirs& dirs,
					  const Momentum& bc,
					  const int& source_or_sink,
					  const bool& include_phases)
  {
    // constexpr std::array<int,2> sign={-1,+1};
    // const int s=sign[source_or_sink]*include_phases;
    
    //compute the main part of the fft
    CRASH("reimplement");
    //fft4d((complex*)out,(complex*)in,dirs,sizeof(spincolor)/sizeof(complex),sign[source_or_sink],0);
    
    // //compute steps
    // Momentum steps;
    // for(int mu=0;mu<NDIM;mu++)
    //   steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    // //add the fractional phase
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	//compute phase exponent
    // 	double arg=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
    // 	//compute the phase
    // 	complex ph={cos(arg),sin(arg)};
	
    // 	//adapt the phase
    // 	for(int id=0;id<NDIRAC;id++)
    // 	  for(int ic=0;ic<NCOL;ic++)
    // 	    safe_complex_prod(out[ivol][id][ic],out[ivol][id][ic],ph);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    
    // set_borders_invalid(out);
  }
  
  void pass_spincolor_from_x_to_mom_space(LxField<spincolor>& out,
					  const LxField<spincolor>& in,
					  const WhichDirs& dirs,
					  const Momentum& bc,
					  const int& source_or_sink,
					  const bool& include_phases)
  {
    
    // constexpr std::array<int,2> sign={+1,-1};
    // const int s=sign[source_or_sink]*include_phases;
    
    // //compute steps
    // Momentum steps;
    // for(int mu=0;mu<NDIM;mu++)
    //   steps[mu]=dirs[mu]*s*bc[mu]*M_PI/glbSize[mu];
    
    // //add the fractional phase
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	//compute phase exponent
    // 	double arg=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  arg+=steps[mu]*glbCoordOfLoclx[ivol][mu];
	
    // 	//compute the phase
    // 	complex ph={cos(arg),sin(arg)};
	
    // 	//adapt the phase
    // 	for(int id=0;id<NDIRAC;id++)
    // 	  for(int ic=0;ic<NCOL;ic++)
    // 	    safe_complex_prod(out[ivol][id][ic],in[ivol][id][ic],ph);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
    
    //compute the main part of the fft
    CRASH("reimplement");
    //fft4d((complex*)out,(complex*)out,dirs,sizeof(spincolor)/sizeof(complex),sign[source_or_sink],1);
  }
}
