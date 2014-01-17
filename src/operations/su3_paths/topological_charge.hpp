#ifndef _TOPOLOGICAL_CHARGE_H
#define _TOPOLOGICAL_CHARGE_H

namespace nissa
{
  void average_topological_charge_eo_conf(double *charge,quad_su3 **eo_conf);
  void average_topological_charge_lx_conf(double *charge,quad_su3 *lx_conf);
  void measure_topology_eo_conf(top_meas_pars_t &pars,quad_su3 **uncooled_conf,int iconf,bool conf_created);
  void measure_topology_lx_conf(top_meas_pars_t &pars,quad_su3 *uncooled_conf,int iconf,bool conf_created,bool presereve_uncooled=true);
  void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf);
  void four_leaves(as2t_su3 *Pmunu,quad_su3 *conf);
  void local_topological_charge(double *charge,quad_su3 *conf);
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in);
  void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in);
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in);
  void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in);
}

#endif
