#ifndef _STOCHASTIC_QQG_VERTEX_H
#define _STOCHASTIC_QQG_VERTEX_H

#include "../types/types.hpp"

void stochastic_x_space_qqg_vertex_source(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false);
void stochastic_x_space_qqg_vertex_source(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false);
void stochastic_x_space_qqg_vertex(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false);
void stochastic_x_space_qqg_vertex(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false);

#endif
