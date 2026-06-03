#ifndef _DIRAC_OPERATOR_TMCLOVQ2_HPP
#define _DIRAC_OPERATOR_TMCLOVQ2_HPP

#include "base/field.hpp"
#include <dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp>

namespace nissa
{
  inline void apply_tmclovQ2(LxField<spincolor>& out,
			     const LxField<quad_su3>& conf,
			     const double& kappa,
			     const std::optional<double>& anis,
			     const LxField<clover_term_t>& Cl,
			     LxField<spincolor>& temp,
			     const double& mu,
			     const LxField<spincolor>& in)
  {
    apply_tmclovQ(temp,conf,kappa,anis,Cl,+mu,in);
    apply_tmclovQ(out,conf,kappa,anis,Cl,-mu,temp);
  }
  
  inline void apply_tmclovQ2_m2(LxField<spincolor>& out,
				const LxField<quad_su3>& conf,
				const double& kappa,
				const std::optional<double>& anis,
				const LxField<clover_term_t>& Cl,
				LxField<spincolor>& temp,
				const double& mu2,
				const LxField<spincolor>& in)
  {
    apply_tmclovQ2(out,conf,kappa,anis,Cl,temp,sqrt(mu2),in);
  }
  
  /// Wraps the application of apply_tmclovQ2_m2 with a specific conf, kappa, Cl and temp
  struct ApplyTmclovQ2M2Functor
  {
    const LxField<quad_su3>& conf;
    
    const double kappa;
    
    const std::optional<double> anis;
    
    const LxField<clover_term_t>& Cl;
    
    LxField<spincolor>& temp;
    
    /// Constructor
    ApplyTmclovQ2M2Functor(const LxField<quad_su3>& conf,
			   const double& kappa,
			   const std::optional<double>& anis,
			   const LxField<clover_term_t>& Cl,
			   LxField<spincolor>& temp) :
      conf(conf),
      kappa(kappa),
      anis(anis),
      Cl(Cl),
      temp(temp)
    {
    }
    
    /// Callable
    void operator()(LxField<spincolor>& out,
		    const double& mass2,
		    const LxField<spincolor>& in)
    {
      apply_tmclovQ2_m2(out,conf,kappa,anis,Cl,temp,mass2,in);
    }
  };
}

#endif
