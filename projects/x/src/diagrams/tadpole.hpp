#ifndef _TADPOLE_H
#define _TADPOLE_H
void compute_tadpole_diagram_in_mom_space(spinspin *q_tad,quark_info qu,gluon_info gl);
void compute_tadpole_diagram_in_x_space(spinspin *q_tad,spin1prop g_prop);
void compute_tadpole_diagram_in_x_space(spinspin *q_tad,quark_info qu,gluon_info gl);
void compute_tadpole_twisted_propagator_in_x_space(spinspin *q_out,quark_info qu,gluon_info gl);
void compute_tadpole_twisted_propagator_in_mom_space(spinspin *q_out,quark_info qu,gluon_info gl);
#endif
