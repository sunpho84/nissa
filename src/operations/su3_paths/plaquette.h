#ifndef _PLAQUETTE_H
#define _PLAQUETTE_H

#include "../../new_types/new_types_definitions.h"

double global_plaquette_eo_conf(quad_su3 **conf);
double global_plaquette_eo_conf_edges(quad_su3 **conf);
double global_plaquette_lx_conf(quad_su3 *conf);
double global_plaquette_variance_lx_conf(quad_su3 *conf);
void compute_point_staples_eo_conf(quad_su3 staple,quad_su3 **eo_conf,int A);
void compute_point_staples_eo_conf_single_dir(su3 staple,quad_su3 **eo_conf,int A,int mu);
void compute_rectangle_staples_eo_conf(quad_su3 **staple,quad_su3 **eo_conf);
void global_plaquette_eo_conf(double *totplaq,quad_su3 **conf);
void global_plaquette_lx_conf(double *totplaq,quad_su3 *conf);
void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu);

#endif
