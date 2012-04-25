#ifndef _GAUGE_FIXING_H
#define _GAUGE_FIXING_H

double compute_landau_or_coulomb_gauge_fixing_quality(quad_su3 *conf,int nmu);
void compute_landau_or_coulomb_delta(su3 g,quad_su3 *conf,int ivol,int nmu);
void compute_landau_or_coulomb_quality_delta(su3 g,quad_su3 *conf,int ivol,int nmu);
void coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision);
void exponentiate(su3 g,su3 a);
void find_landau_or_coulomb_gauge_fixing_matr(su3 *fixm,quad_su3 *conf,double required_precision,int nmu);
void find_local_landau_or_coulomb_gauge_fixing_transformation(su3 g,quad_su3 *conf,int ivol,int nmu);
void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u);
void gauge_transform_conf(quad_su3 *uout,su3 *g,quad_su3 *uin);
void landau_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision);
void landau_or_coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision,int nmu);
void local_gauge_transform(quad_su3 *conf,su3 g,int ivol);
void overrelax(su3 out,su3 in,double w);

#endif
