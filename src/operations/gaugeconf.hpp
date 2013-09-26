#ifndef _gaugeconf_hpp
#define _gaugeconf_hpp

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  inline int check_add_square_staple(int *isquare_staples_to_ask,int &nsquare_staple_to_ask,int ivol,int dir,int verse,int iter);
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
  void ac_rotate_vector(void *out,void *in,int axis,int bps);
  void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges);
  void cool_conf(quad_su3 **eo_conf,int over_flag,double over_exp);
  void generate_cold_eo_conf(quad_su3 **conf);
  void generate_hot_eo_conf(quad_su3 **conf);
  void heatbath_conf(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars);
  void heatbath_or_overrelax_conf(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars,int heat_over);
  void heatbath_or_overrelax_conf_Wilson_action(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars,int heat_over);
  void heatbath_or_overrelax_conf_tlSym_action(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars,int heat_over);
  void overrelax_conf(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars);
  void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void unitarity_check_lx_conf(unitarity_check_result_t &result,quad_su3 *conf);
}

#endif
