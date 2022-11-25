#ifndef _TWISTED_DIRAC_EOPREC_OPERATOR_HPP
#define _TWISTED_DIRAC_EOPREC_OPERATOR_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
#include <free_theory/free_theory_types.hpp>
#include <new_types/spin.hpp>

namespace nissa
{
  void tmDkern_eoprec_eos(OddField<spin>& out,
			  EvnField<spin>& temp,
			  const tm_quark_info& qu,
			  const OddField<spin>& in);
  
  void tmDkern_eoprec_square_eos(OddField<spin>& out,
				 OddField<spin>& temp1,
				 EvnField<spin> &temp2,
				 const tm_quark_info& qu,
				 const OddField<spin>& in);
  
  /// Inverse
  template <typename O,
	    typename I>
  void inv_tmDee_or_oo_eos(O&& out,
			   const tm_quark_info& qu,
			   const I& in)
  {
    const double a=1/(2*qu.kappa),b=qu.mass,nrm=1/(a*a+b*b);
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	const complex z={+a*nrm,-b*nrm};
	for(int id=0;id<NDIRAC/2;id++)
	  unsafe_complex_prod(out[X][id],in[X][id],z);
	for(int id=NDIRAC/2;id<NDIRAC;id++)
	  unsafe_complex_conj2_prod(out[X][id],in[X][id],z);
      }
    NISSA_PARALLEL_LOOP_END;
    
    out.invalidateHalo();
  }
  
  //apply even-odd or odd-even part of tmD, multiplied by -2
  template <typename O,
	    typename I>
  void tmn2Deo_or_tmn2Doe_eos(O& out,
			      const I& in,
			      const momentum_t& bc)
  {
    in.updateHalo();
    
    std::array<complex,NDIM> phases;
    for(int mu=0;mu<NDIM;mu++)
      {
	phases[mu][RE]=cos(M_PI*bc[mu]);
	phases[mu][IM]=sin(M_PI*bc[mu]);
      }
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	spin_put_to_zero(out[X]);
	
	for(int mu=0;mu<NDIM;mu++)
	  {
	    complex temp_c0,temp_c1;
	    
	    //Forward
	    const int Xup=I::locNeighup(X,mu);
	    switch(mu)
	    {
	    case 0:
	      complex_summ(temp_c0,in[Xup][0],in[Xup][2]);
	      complex_summ(temp_c1,in[Xup][1],in[Xup][3]);
	      break;
	    case 1:
	      complex_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	      complex_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	      break;
	    case 2:
	      complex_summ(temp_c0,in[Xup][0],in[Xup][3]);
	      complex_subt(temp_c1,in[Xup][1],in[Xup][2]);
	      break;
	    case 3:
	      complex_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	      complex_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	      break;
	    }
	    
	    if(I::glbCoord(Xup,mu)==0)
	      {
		safe_complex_prod(temp_c0,temp_c0,phases[mu]);
		safe_complex_prod(temp_c1,temp_c1,phases[mu]);
	      }
	    
	    complex_summassign(out[X][0],temp_c0);
	    complex_summassign(out[X][1],temp_c1);
	    
	    switch(mu)
	      {
	      case 0:
		complex_summassign(out[X][2],out[X][0]);
		complex_summassign(out[X][3],out[X][1]);
		break;
	      case 1:
		complex_isubtassign(out[X][2],temp_c1);
		complex_isubtassign(out[X][3],temp_c0);
		break;
	      case 2:
		complex_subtassign(out[X][2],temp_c1);
		complex_summassign(out[X][3],temp_c0);
		break;
	      case 3:
		complex_isubtassign(out[X][2],temp_c0);
		complex_isummassign(out[X][3],temp_c1);
		break;
	      }
	    
	    //Backward
	    const int Xdw=I::locNeighdw(X,mu);
	    switch(mu)
	    {
	    case 0:
	      complex_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	      complex_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	      break;
	    case 1:
	      complex_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	      complex_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	      break;
	    case 2:
	      complex_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	      complex_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	      break;
	    case 3:
	      complex_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	      complex_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	      break;
	    }
	    
	    if(O::glbCoord(X,mu)==0)
	      {
		safe_complex_conj2_prod(temp_c0,temp_c0,phases[mu]);
		safe_complex_conj2_prod(temp_c1,temp_c1,phases[mu]);
	      }
	    
	    complex_summassign(out[X][0],temp_c0);
	    complex_summassign(out[X][1],temp_c1);
	    
	    switch(mu)
	      {
	      case 0:
		complex_subtassign(out[X][2],temp_c0);
		complex_subtassign(out[X][3],temp_c1);
		break;
	      case 1:
		complex_isummassign(out[X][2],temp_c1);
		complex_isummassign(out[X][3],temp_c0);
		break;
	      case 2:
		complex_summassign(out[X][2],temp_c1);
		complex_subtassign(out[X][3],temp_c0);
		break;
	      case 3:
		complex_isummassign(out[X][2],temp_c0);
		complex_isubtassign(out[X][3],temp_c1);
		break;
	      }
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
}

#endif
