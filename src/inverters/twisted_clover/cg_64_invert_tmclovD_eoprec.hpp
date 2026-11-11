#ifndef _CG_64_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_64_INVERT_TMCLOVD_EOPREC_HPP

#include <optional>

#include <dirac_operators/tmQ/dirac_operator_tmQ.hpp>
#include <dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp>
#include <inverters/templates/cg_invert_template_threaded.hpp>

namespace nissa
{
  inline void inv_tmclovDkern_eoprec_square_eos_cg_64_portable(OddField<spincolor>& sol,
							       std::optional<OddField<spincolor>> guess,
							       const EoField<quad_su3>& conf,
							       const double& kappa,
							       const AnisDopPars& anisDopPars,
							       const OddField<clover_term_t>& Cl_odd,
							       const EvnField<inv_clover_term_t>& invCl_evn,
							       const double& mu,
							       const int& niter,
							       const double residue,
							       const OddField<spincolor>& source)
  {
    std::function<void(OddField<spincolor>&,
      const OddField<spincolor>&)> f=
      [temp1=OddField<spincolor>("temp1",WITH_HALO),
       temp2=EvnField<spincolor>("temp2",WITH_HALO),
       &conf,
       &kappa,
       &anisDopPars,
       &Cl_odd,
       &invCl_evn,
       &mu](OddField<spincolor>& out,
	    const OddField<spincolor>& in) mutable
      {
	tmclovDkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,anisDopPars,Cl_odd,invCl_evn,mu,in);
      };
    
    cg_invert(sol,
	      guess,
	      f,
	      niter,
	      residue,
	      source);
  }
  
  /// Run native CG after external solver selection in the full clover inverter
  inline void inv_tmclovDkern_eoprec_square_eos_cg_64(OddField<spincolor>& sol,
						      std::optional<OddField<spincolor>> guess,
						      const EoField<quad_su3>& conf,
						      const double& kappa,
						      const AnisDopPars& anisDopPars,
						      const OddField<clover_term_t>& Cl_odd,
						      const EvnField<inv_clover_term_t>& invCl_evn,
						      const double& mu,
						      const int& niter,
						      const double residue,
						      const OddField<spincolor>& source)
  {
    inv_tmclovDkern_eoprec_square_eos_cg_64_portable(sol,guess,conf,kappa,anisDopPars,Cl_odd,invCl_evn,mu,niter,residue,source);
  }
}

#endif
