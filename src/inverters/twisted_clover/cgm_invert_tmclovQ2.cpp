#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp"
#include "inverters/templates/cgm_invert_template_threaded.hpp"

namespace nissa
{
  void inv_tmclovQ2_cgm(std::vector<LxField<spincolor>>& sol,
			const LxField<quad_su3>& conf,
			const double& kappa,
			const LxField<clover_term_t>& Cl,
			const std::vector<double>& m,
			const int& niter_max,
			const std::vector<double>& req_res,
			const LxField<spincolor>& source)
  {
    std::vector<double> m2(m);
    for(auto& m : m2)
      m*=m;
    
    cgm_invert(sol,
	       m2,
	       [temp=LxField<spincolor>("temp",WITH_HALO),
		&conf,
		&kappa,
		&Cl]
	       (LxField<spincolor>& out,
		const double& mass2,
		const LxField<spincolor>& in) mutable
	       {
		 apply_tmclovQ2_m2(out,conf,kappa,Cl,temp,mass2,in);
	       },
	       niter_max,
	       req_res,
	       source);
  }
}
