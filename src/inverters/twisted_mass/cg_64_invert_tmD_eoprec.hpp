#ifndef _CG_64_INVERT_TMDEOIMPR_HPP
#define _CG_64_INVERT_TMDEOIMPR_HPP

#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp>
#include <dirac_operators/tmQ/dirac_operator_tmQ.hpp>
#include <inverters/templates/cg_invert_template_threaded.hpp>

namespace nissa
{
  inline void inv_tmDkern_eoprec_square_eos_cg_64(OddField<spincolor>& sol,
						  std::optional<OddField<spincolor>> guess,
						  const EoField<quad_su3>& conf,
						  const double& kappa,
						  const AnisDopPars& anisDopPars,
						  const double& mu,
						  const int& niter,
						  const double residue,
						  const OddField<spincolor>& source)
  {
    std::function<void(OddField<spincolor>& out,
		    const OddField<spincolor>& in)> f=
      [temp1=OddField<spincolor>("temp1",WITH_HALO),
       temp2=EvnField<spincolor>("temp2",WITH_HALO),
       &conf,
       &kappa,
       &anisDopPars,
       &mu](OddField<spincolor>& out,
	    const OddField<spincolor>& in) mutable
      {
	tmDkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,anisDopPars,mu,in);
      };
    
    cg_invert(sol,
	      guess,
	      f,
	      niter,
	      residue,
	      source);
  }
}

#endif
