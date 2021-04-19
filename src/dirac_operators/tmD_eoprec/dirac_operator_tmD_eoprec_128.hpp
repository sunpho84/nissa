#ifndef _DIRAC_OPERATOR_TMDEOIMPR_128_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_128_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/float_128.hpp"

namespace nissa
{
  void tmn2Deo_or_tmn2Doe_eos_128(spincolor_128 *out,eo_ptr<quad_su3> conf,const Parity& eooe,spincolor_128 *in);
  //wrappers
  inline void tmn2Doe_eos_128(spincolor_128 *out,eo_ptr<quad_su3> conf,spincolor_128 *in){tmn2Deo_or_tmn2Doe_eos_128(out,conf,1,in);}//nasty
  inline void tmn2Deo_eos_128(spincolor_128 *out,eo_ptr<quad_su3> conf,spincolor_128 *in){tmn2Deo_or_tmn2Doe_eos_128(out,conf,0,in);}
  
  void inv_tmDee_or_oo_eos_128(spincolor_128 *out,double kappa,double mu,spincolor_128 *in);
  void tmDee_or_oo_eos_128(spincolor_128 *out,double kappa,double mu,spincolor_128 *in);
  void tmDkern_eoprec_eos_128(spincolor_128 *out,spincolor *temp,eo_ptr<quad_su3> conf,double kappa,double mu,spincolor_128 *in);
  void tmDkern_eoprec_square_eos_128(spincolor_128 *out,spincolor_128 *temp1,spincolor_128 *temp2,eo_ptr<quad_su3> conf,double kappa,double mu,spincolor_128 *in);
  void tmDkern_eoprec_eos_put_together_and_include_gamma5_128(spincolor_128 *out,spincolor_128 *temp);
}

#endif
