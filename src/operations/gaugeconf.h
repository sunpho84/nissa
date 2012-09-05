#ifndef _GAUGECONF_H
#define _GAUGECONF_H

void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
void ac_rotate_vector(void *out,void *in,int axis,int bps);
void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges);
void generate_cold_eo_conf(quad_su3 **conf);
void generate_hot_eo_conf(quad_su3 **conf);
void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
void cool_conf(quad_su3 **eo_conf,int over_flag,double over_exp);

#endif
