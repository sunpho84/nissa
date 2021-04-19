#ifndef _TWISTED_DIRAC_EOPREC_OPERATOR_HPP
#define _TWISTED_DIRAC_EOPREC_OPERATOR_HPP

#include <new_types/spin.hpp>
#include <free_theory/free_theory_types.hpp>

namespace nissa
{
  /// Apply even-odd or odd-even part of tmD, multiplied by -2
  void tmn2Deo_or_tmn2Doe_eos(spin *out,int eooe,spin *in,const Momentum& bc);
  
  void tmn2Doe_eos(spin *out,spin *in,const Momentum& bc);
  void tmn2Deo_eos(spin *out,spin *in,const Momentum& bc);
  void tmDee_or_oo_eos(spin *out,const tm_quark_info& qu,spin *in);
  void inv_tmDee_or_oo_eos(spin *out,const tm_quark_info& qu,spin *in);
  void tmDkern_eoprec_eos(spin *out,spin *temp,const tm_quark_info& qu,spin *in);
  void tmDkern_eoprec_square_eos(spin *out,spin *temp1,spin *temp2,const tm_quark_info& qu,spin *in);
}

#endif
