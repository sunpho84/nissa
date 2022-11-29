#ifndef _DIRAC_OPERATOR_TMD_EOIMPR_HPP
#define _DIRAC_OPERATOR_TMD_EOIMPR_HPP

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //improve
  
  // void inv_tmDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  // void tmDee_or_oo_eos(spincolor *out,double kappa,double mu,spincolor *in);
  // void tmDkern_eoprec_eos(spincolor *out,spincolor *temp,eo_ptr<quad_su3>  conf,double kappa,double mu,spincolor *in);
  void tmDkern_eoprec_square_eos(OddField<spincolor>& out,
				 OddField<spincolor>& temp1,
				 EvnField<spincolor> &temp2,
				 EoField<quad_su3>& conf,
				 const double& kappa,
				 const double& mu,
				 const OddField<spincolor>& in);
  // void tmn2Deo_or_tmn2Doe_eos(spincolor *out,eo_ptr<quad_su3> conf,int eooe,spincolor *in);
  // void tmDkern_eoprec_eos_put_together_and_include_gamma5(spincolor *out,spincolor *temp);
}

#endif
