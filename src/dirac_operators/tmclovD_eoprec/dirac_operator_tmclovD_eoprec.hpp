#ifndef _DIRAC_OPERATOR_TMCLOVD_EOPREC_HPP
#define _DIRAC_OPERATOR_TMCLOVD_EOPREC_HPP

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void tmclovDee_or_oo_eos(spincolor *out,double kappa,clover_term_t *Cl,bool dag,double mu,spincolor *in);
  void inv_tmclovDee_or_oo_eos(spincolor *out,inv_clover_term_t *invCl,bool dag,spincolor *in);
  void tmclovDkern_eoprec_eos(spincolor *out,spincolor *temp,eo_ptr<quad_su3> conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,bool dag,double mu,spincolor *in);
  
  void tmclovDkern_eoprec_square_eos(OddField<spincolor>& out,
				     OddField<spincolor>& temp1,
				     EvenOrOddField<spincolor>& temp2,
				     const EoField<quad_su3>& conf,
				     const double& kappa,
				     const OddField<clover_term_t>& Cl_odd,
				     const EvnField<inv_clover_term_t>& invCl_evn,
				     const double& mu,
				     const OddField<spincolor>& in);
}

#endif
