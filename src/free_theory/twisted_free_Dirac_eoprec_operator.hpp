#ifndef _TWISTED_DIRAC_EOPREC_OPERATOR_HPP
#define _TWISTED_DIRAC_EOPREC_OPERATOR_HPP

#include "new_types/spin.hpp"
#include "free_theory_types.hpp"

namespace nissa
{
  void tmn2Deo_or_tmn2Doe_eos(spin *out,int eooe,spin *in,const momentum_t& bc);
  void tmn2Doe_eos(spin *out,spin *in,const momentum_t& bc);
  void tmn2Deo_eos(spin *out,spin *in,const momentum_t& bc);
  void tmDee_or_oo_eos(spin *out,tm_quark_info qu,spin *in);
  void inv_tmDee_or_oo_eos(spin *out,tm_quark_info qu,spin *in);
  void tmDkern_eoprec_eos(spin *out,spin *temp,tm_quark_info qu,spin *in);
  void tmDkern_eoprec_square_eos(spin *out,spin *temp1,spin *temp2,tm_quark_info qu,spin *in);
}

#endif
