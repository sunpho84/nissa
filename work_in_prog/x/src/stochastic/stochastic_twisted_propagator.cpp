#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/vectors.h"
#include "../../../../src/operations/fft.h"

#include "../routines/fourier.h"
#include "../types/types.h"
#include "../vertex/x_space_stochastic_qqg_vertex.h"

void generate_stochastic_A_twisted_propagator_by_inv(spinspin *q_A,spinspin *q_prop,quark_info qu,spin1field *g_A,gluon_info gl)
{stochastic_x_space_qqg_vertex(q_A,q_prop,qu,g_A,gl);}

void generate_stochastic_A_dag_twisted_propagator_by_inv(spinspin *q_A,spinspin *q_prop,quark_info qu,spin1field *g_A,gluon_info gl)
{stochastic_x_space_qqg_vertex(q_A,q_prop,qu,g_A,gl,true);}

void generate_stochastic_A_B_dag_twisted_propagator_by_inv(spinspin *q_AB,spinspin *q_prop,quark_info qu,spin1field *g_A,spin1field *g_B,gluon_info gl)
{
  stochastic_x_space_qqg_vertex(q_AB,q_prop,qu,g_A,gl);
  stochastic_x_space_qqg_vertex(q_AB,q_AB,qu,g_B,gl,true);
}

void generate_stochastic_A_B_dag_twisted_propagator_source_by_inv(spinspin *q_AB,spinspin *q_prop,quark_info qu,spin1field *g_A,spin1field *g_B,gluon_info gl)
{
  stochastic_x_space_qqg_vertex(q_AB,q_prop,qu,g_A,gl);
  stochastic_x_space_qqg_vertex_source(q_AB,q_AB,qu,g_B,gl,true);
}
