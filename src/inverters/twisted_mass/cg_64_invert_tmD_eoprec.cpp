#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp>
#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.hpp>
#include <inverters/templates/cg_invert_template_threaded.hpp>

namespace nissa
{
  void inv_tmDkern_eoprec_square_eos_cg_64(OddField<spincolor>& sol,
					   std::optional<OddField<spincolor>> guess,
					   const EoField<quad_su3>& conf,
					   const double& kappa,
					   const double& mu,
					   const int& niter,
					   const double residue,
					   const OddField<spincolor>& source)
  {
    OddField<spin> temp1("temp1",WITH_HALO);
    EvnField<spin> temp2("temp2",WITH_HALO);
    
    cg_invert(sol,
	      guess,
	      [temp1=OddField<spincolor>("temp1",WITH_HALO),
	       temp2=EvnField<spincolor>("temp2",WITH_HALO),
	       &conf,
	       &kappa,
	       &mu](OddField<spincolor>& out,
		    const OddField<spincolor>& in) mutable
	      {
		tmDkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,mu,in);
	      },
	      niter,
	      residue,
	      source);
  }
}
