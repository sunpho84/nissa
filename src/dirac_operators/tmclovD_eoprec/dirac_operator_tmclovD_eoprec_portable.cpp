#pragma once

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "threads/threads.hpp"

#include <dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp>
#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.hpp>

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //implement ee or oo part of Dirac operator, equation(3)
  /// Inverse
  template <typename O,
	    typename C,
	    typename I>
  void tmclovDee_or_oo_eos(O&& out,
			   const double& kappa,
			   const C& Cl,
			   const bool& dag,
			   double mu,
			   const I& in)
  {
    if(dag) mu=-mu;
    
    if(in==out) crash("in==out!");
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	apply_point_twisted_clover_term_to_halfspincolor(out[X],+mu,kappa,Cl[X],in[X],0);
	apply_point_twisted_clover_term_to_halfspincolor(out[X],-mu,kappa,Cl[X],in[X],NDIRAC/2);
      }
    NISSA_PARALLEL_LOOP_END;
    
    out.invalidateHalo();
  }
  
  //implement Koo defined in equation (7)
  void tmclovDkern_eoprec_eos(OddField<spincolor>& out,
			      EvnField<spincolor>& tmp,
			      const EoField<quad_su3>& conf,
			      const double& kappa,
			      const OddField<clover_term_t>& Cl_odd,
			      const EvnField<inv_clover_term_t>& invCl_evn,
			      const bool& dag,
			      const double& mu,
			      const OddField<spincolor>& in)
  {
    /// Improve
    EvnField<spincolor>& outAsEvnTmp=out.castSitesCoverage<EVEN_SITES>();
    // EvnField<spincolor>& tmpEvn=extTmp.castSitesCoverage<EVEN_SITES>();
    // OddField<spincolor>& tmpOdd=extTmp.castSitesCoverage<ODD_SITES>();
    tmn2Deo_or_tmn2Doe_eos(outAsEvnTmp,conf,in);
    inv_tmclovDee_or_oo_eos(tmp.castSitesCoverage<EVEN_SITES>(),invCl_evn,dag,outAsEvnTmp);
    tmn2Deo_or_tmn2Doe_eos(out,conf,tmp.castSitesCoverage<EVEN_SITES>());
    
    tmclovDee_or_oo_eos(tmp.castSitesCoverage<ODD_SITES>(),kappa,Cl_odd,dag,mu,in);
    
    tmDkern_eoprec_eos_put_together_and_include_gamma5(out,tmp.castSitesCoverage<ODD_SITES>());
  }
  
  //square of Koo
  void tmclovDkern_eoprec_square_eos(OddField<spincolor>& out,
				     OddField<spincolor>& temp1,
				     EvnField<spincolor>& temp2,
				     const EoField<quad_su3>& conf,
				     const double& kappa,
				     const OddField<clover_term_t>& Cl_odd,
				     const EvnField<inv_clover_term_t>& invCl_evn,
				     const double& mu,
				     const OddField<spincolor>& in)
  {
    tmclovDkern_eoprec_eos(temp1,temp2,conf,kappa,Cl_odd,invCl_evn,true,  mu,in   );
    tmclovDkern_eoprec_eos(out,  temp2,conf,kappa,Cl_odd,invCl_evn,false, mu,temp1);
  }
}
