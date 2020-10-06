#include "../../../../src/nissa.hpp"
using namespace std;
using namespace nissa;

#include "../types/types.hpp"
#include "../vertex/x_space_stochastic_qqg_vertex.hpp"
#include "../propagators/twisted_propagator.hpp"

void generate_stochastic_A_twisted_propagator(spinspin *q_A,spinspin *q_prop,quark_info qu,spin1field *g_A,gluon_info gl)
{stochastic_x_space_qqg_vertex(q_A,q_prop,qu,g_A,gl);}

void generate_stochastic_A_dag_twisted_propagator(spinspin *q_A,spinspin *q_prop,quark_info qu,spin1field *g_A,gluon_info gl)
{stochastic_x_space_qqg_vertex(q_A,q_prop,qu,g_A,gl,true);}


void generate_stochastic_A_B_dag_twisted_propagator_source(spinspin *q_AB,spinspin *q_prop,quark_info qu,spin1field *g_A,spin1field *g_B,gluon_info gl)
{
  spinspin *q_B=nissa_malloc("q_B",loc_vol+bord_vol,spinspin);
  
  stochastic_x_space_qqg_vertex(q_B,q_prop,qu,g_B,gl,true);
  stochastic_x_space_qqg_vertex_source(q_AB,q_B,qu,g_A,gl);
  
  nissa_free(q_B);
}

void generate_stochastic_A_B_dag_twisted_propagator(spinspin *q_AB,spinspin *q_prop,quark_info qu,spin1field *g_A,spin1field *g_B,gluon_info gl)
{
  generate_stochastic_A_B_dag_twisted_propagator_source(q_AB,q_prop,qu,g_A,g_B,gl);
  multiply_from_left_by_x_space_twisted_propagator_by_fft(q_AB,q_AB,qu);
}
