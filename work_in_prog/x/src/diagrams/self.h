#ifndef _SELF_H
#define _SELF_H

#include "../../../../src/new_types/new_types_definitions.h"

#include "../types/types.h"

void compute_self_energy_twisted_propagator_in_mom_space(spinspin *q_out,quark_info qu,gluon_info gl);
void compute_self_energy_twisted_propagator_in_mom_space(spinspin *q_out,spinspin *q_prop,quark_info qu,spin1prop *g_prop,gluon_info gl);
void compute_self_energy_twisted_propagator_in_x_space(spinspin *q_out,quark_info qu,gluon_info gl);
void compute_self_energy_twisted_propagator_in_x_space(spinspin *q_out,spinspin *q_prop,quark_info qu,spin1prop *g_prop,gluon_info gl);
void summ_the_contribution_of_self_energy_twisted_propagator_in_x_space(spinspin *q_out,spinspin *osi,spinspin *q,spin1prop *g,int nu,int mu,spinspin *oso,double sign);

#endif
