#ifndef _GAUGECONF_HPP
#define _GAUGECONF_HPP

#include "operations/su3_paths/gauge_sweeper.hpp"
#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  inline int check_add_square_staple(int *isquare_staples_to_ask,int &nsquare_staple_to_ask,int ivol,int dir,int verse,int iter);
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
  void ac_rotate_vector(void *out,void *in,int axis,size_t bps);
  void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges);
  void generate_cold_eo_conf(quad_su3 **conf);
  void generate_hot_eo_conf(quad_su3 **conf);
  void generate_cold_lx_conf(quad_su3 *conf);
  void generate_hot_lx_conf(quad_su3 *conf);
  void heatbath_lx_conf(gauge_sweeper_t *sweeper,quad_su3 *conf,double beta,int nhits);
  void overrelax_lx_conf(gauge_sweeper_t *sweeper,quad_su3 *conf,int nhits);
  void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void unitarity_check_lx_conf(unitarity_check_result_t &result,quad_su3 *conf);
  void unitarize_lx_conf_orthonormalizing(quad_su3 *conf);
  void unitarize_lx_conf_maximal_trace_projecting(quad_su3 *conf);
  void unitarize_eo_conf_maximal_trace_projecting(quad_su3 **conf);
}

#endif
