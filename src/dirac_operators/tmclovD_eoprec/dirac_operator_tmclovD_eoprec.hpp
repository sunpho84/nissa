#ifndef _DIRAC_OPERATOR_TMCLOVD_EOPREC_HPP
#define _DIRAC_OPERATOR_TMCLOVD_EOPREC_HPP

#include "base/field.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  void tmclovDee_or_oo_eos(spincolor *out,double kappa,clover_term_t *Cl,bool dag,double mu,spincolor *in);
  
  void tmclovDkern_eoprec_eos(OddField<spincolor>& out,
			      EvnField<spincolor>& tmp,
			      const EoField<quad_su3>& conf,
			      const double& kappa,
			      const OddField<clover_term_t>& Cl_odd,
			      const EvnField<inv_clover_term_t>& invCl_evn,
			      const bool& dag,
			      const double& mu,
			      const OddField<spincolor>& in);
  
  //inverse
  template <typename O,
	    typename C,
	    typename I>
  void inv_tmclovDee_or_oo_eos(O&& out,
			       const C& invCl,
			       const bool& dag,
			       const I& in)
  {
    //if dagger, swaps the sign of mu, which means taking the hermitian of the inverse
    int high=0,low=1;
    if(dag) std::swap(low,high);
    
    PAR(0,locVolh,
	CAPTURE(high,low,
		TO_WRITE(out),
		TO_READ(in),
		TO_READ(invCl)),X,
	{
	  unsafe_halfspincolor_halfspincolor_times_halfspincolor(out[X],invCl[X][high],in[X],2*high);
	  unsafe_halfspincolor_halfspincolor_dag_times_halfspincolor(out[X],invCl[X][low],in[X],2*low);
	});
  }
  
  void tmclovDkern_eoprec_square_eos(OddField<spincolor>& out,
				     OddField<spincolor>& temp1,
				     EvnField<spincolor>& temp2,
				     const EoField<quad_su3>& conf,
				     const double& kappa,
				     const OddField<clover_term_t>& Cl_odd,
				     const EvnField<inv_clover_term_t>& invCl_evn,
				     const double& mu,
				     const OddField<spincolor>& in);
}

#endif
