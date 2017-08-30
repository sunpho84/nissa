#ifndef _GAUGE_FIXING_HPP
#define _GAUGE_FIXING_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  double compute_Landau_or_Coulomb_gauge_fixing_quality(quad_su3 *conf,int nmu);
  void compute_Landau_or_Coulomb_delta(su3 g,quad_su3 *conf,int ivol,int nmu);
  void compute_Landau_or_Coulomb_quality_delta(su3 g,quad_su3 *conf,int ivol,int nmu);
  void Coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision);
  void exponentiate(su3 g,su3 a);
  void find_Landau_or_Coulomb_gauge_fixing_matr(su3 *fixm,quad_su3 *conf,double required_precision,int nmu);
  void find_local_Landau_or_Coulomb_gauge_fixing_transformation(su3 g,quad_su3 *conf,int ivol,int nmu);
  void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u);
  void gauge_transform_conf(quad_su3 *uout,su3 *g,const quad_su3 *uin);
  void gauge_transform_conf(quad_su3 **uout,su3 **g,quad_su3 **uin);
  void gauge_transform_color(color **out,su3 **g,color **in);
  void Landau_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision);
  void Landau_or_Coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision,int nmu);
  void local_gauge_transform(quad_su3 *conf,su3 g,int ivol);
  void perform_random_gauge_transform(quad_su3 *conf_out,quad_su3 *conf_in);
  void perform_random_gauge_transform(quad_su3 **conf_out,quad_su3 **conf_in);
}

#endif
