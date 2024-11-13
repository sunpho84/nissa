#ifndef _CG_64_INVERT_TMCLOVD_EOPREC_HPP
#define _CG_64_INVERT_TMCLOVD_EOPREC_HPP

#include <optional>

#include <dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp>
#include <inverters/templates/cg_invert_template_threaded.hpp>

namespace nissa
{
  inline void inv_tmclovDkern_eoprec_square_eos_cg_64_portable(OddField<spincolor>& sol,
							       std::optional<OddField<spincolor>> guess,
							       const EoField<quad_su3>& conf,
							       const double& kappa,
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
       &Cl_odd,
       &invCl_evn,
       &mu](OddField<spincolor>& out,
	    const OddField<spincolor>& in) mutable
      {
	tmclovDkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,Cl_odd,invCl_evn,mu,in);
      };
    
    cg_invert(sol,
	      guess,
	      f,
	      niter,
	      residue,
	      source);
  }
  
  //wrapper
  inline void inv_tmclovDkern_eoprec_square_eos_cg_64(OddField<spincolor>& sol,
						      std::optional<OddField<spincolor>> guess,
						      const EoField<quad_su3>& conf,
						      const double& kappa,
						      const OddField<clover_term_t>& Cl_odd,
						      const EvnField<inv_clover_term_t>& invCl_evn,
						      const double& mu,
						      const int& niter,
						      const double residue,
						      const OddField<spincolor>& source)
  {
#if defined USE_DDALPHAAMG
    crash("reimplement");//move upward, reorder everything
    
    // if(use_DD and fabs(mu)<=DD::max_mass)
    // 	{
    // 	  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol,quad_su3);
    // 	  paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    // 	  spincolor *tmp_in=nissa_malloc("tmp_in",loc_vol,spincolor);
    // 	  spincolor *tmp_out=nissa_malloc("tmp_out",loc_vol,spincolor);
    // 	  vector_reset(tmp_in);
    // 	  double_vector_copy((double*)tmp_in,(double*)source,loc_volh*sizeof(spincolor)/sizeof(double));
    // 	  DD::solve(tmp_out,lx_conf,kappa,cSW,mu,residue,tmp_in,true);
    // 	  nissa_free(lx_conf);
    // 	  inv_tmclovDkern_eoprec_square_eos_cg_64_portable(sol,guess,eo_conf,kappa,Cl_odd,invCl_evn,mu,niter,residue,source);
    // 	  master_printf("%lg %lg\n",tmp_out[0][0][0][0],sol[0][0][0][0]);
    // 	  master_printf("%lg %lg\n",tmp_out[0][0][0][1],sol[0][0][0][1]);
    // 	  nissa_free(tmp_out);
    // 	  nissa_free(tmp_in);
    // 	}
    // else
    master_printf("DDalpha not yet working, probably expecting a different layout\n");
#endif
    inv_tmclovDkern_eoprec_square_eos_cg_64_portable(sol,guess,conf,kappa,Cl_odd,invCl_evn,mu,niter,residue,source);
  }
}

#endif
