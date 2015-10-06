#ifndef _TOPOLOGICAL_CHARGE_HPP
#define _TOPOLOGICAL_CHARGE_HPP

namespace nissa
{
  void build_chromo_therm_from_anti_symmetric_four_leaves(quad_su3 out,as2t_su3 in);
  void measure_topology_eo_conf(top_meas_pars_t &pars,quad_su3 **uncooled_conf,int iconf,bool conf_created);
  void measure_topology_lx_conf(top_meas_pars_t &pars,quad_su3 *uncooled_conf,int iconf,bool conf_created,bool presereve_uncooled=false);
  double topodynamical_potential(double Q,topotential_pars_t &pars);
  void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf);
  void opt_Pmunu_term(quad_su3 *C,quad_su3 *conf);
  void draw_topodynamical_potential(topotential_pars_t &pars);
  void four_leaves(as2t_su3 *Pmunu,quad_su3 *conf);
  void local_topological_charge(double *charge,quad_su3 *conf);
  void topological_staples(quad_su3 *staples,quad_su3 *conf);
  void total_topological_charge_eo_conf(double *charge,quad_su3 **eo_conf);
  void total_topological_charge_lx_conf(double *charge,quad_su3 *lx_conf);
  void unsafe_apply_opt_chromo_operator_to_spincolor(spincolor *out,opt_as2t_su3 *Cl,spincolor *in);
  void unsafe_apply_opt_chromo_operator_to_spincolor_128(spincolor_128 *out,opt_as2t_su3 *Cl,spincolor_128 *in);
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in);
  void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in);
  void unsafe_apply_chromo_operator_to_spincolor_128(spincolor_128 *out,as2t_su3 *Pmunu,spincolor_128 *in);
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in);
  void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in);
}

#endif
