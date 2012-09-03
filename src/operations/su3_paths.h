#ifndef _SU3_PATHS_H
#define _SU3_PATHS_H

double average_real_part_of_trace_of_rectangle_path(quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u);
double global_plaquette_eo_conf(quad_su3 **conf);
void global_plaquette_eo_conf(double *plaq,quad_su3 **conf);
double global_plaquette_lx_conf(quad_su3 *conf);
double global_plaquette_variance_lx_conf(quad_su3 *conf);
void Pline(su3 *Pline,quad_su3 *conf);
void Pline_backward(su3 *Pline, quad_su3 *conf);
void Pline_forward(su3 *Pline, quad_su3 *conf);
void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf);
void average_trace_of_rectangle_path(complex tra,quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u);
void compute_point_staples_eo_conf(quad_su3 staple,quad_su3 **eo_conf,int A);
void compute_point_staples_eo_conf_single_dir(quad_su3 staple,quad_su3 **eo_conf,int A,int mu);
void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu);
void su3_vec_single_shift(su3 *u,int mu,int sign);
void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in);
void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in);
void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in);
void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in);
void average_polyakov_loop(complex tra,quad_su3 *eo_conf,int mu);
void average_polyakov_loop_of_eos_conf(complex tra,quad_su3 **eo_conf,int mu);
double average_topological_charge(quad_su3 *conf);
double average_topological_charge(quad_su3 **conf);
void local_topological_charge(double *charge,quad_su3 *conf);

#endif
