#ifndef _twisted_Dirac_eoprec_operator_H
#define _twisted_Dirac_eoprec_operator_H

#include "../../../../src/new_types/new_types_definitions.h"

#include "../types/types.h"

void tmn2Deo_or_tmn2Doe_eos(spin *out,int eooe,spin *in,momentum_t bc);
void tmn2Doe_eos(spin *out,spin *in,momentum_t bc);
void tmn2Deo_eos(spin *out,spin *in,momentum_t bc);
void tmDee_or_oo_eos(spin *out,quark_info qu,spin *in);
void inv_tmDee_or_oo_eos(spin *out,quark_info qu,spin *in);
void tmDkern_eoprec_eos(spin *out,spin *temp,quark_info qu,spin *in);
void tmDkern_eoprec_square_eos(spin *out,spin *temp1,spin *temp2,quark_info qu,spin *in);

#endif
