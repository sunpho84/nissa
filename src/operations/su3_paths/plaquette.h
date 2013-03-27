#ifndef _PLAQUETTE_H
#define _PLAQUETTE_H

#include "../../new_types/new_types_definitions.h"

double global_plaquette_lx_conf(quad_su3 *conf);
double global_plaquette_eo_conf(quad_su3 **conf);
void global_plaquette_eo_conf(double *totplaq,quad_su3 **conf);
void global_plaquette_lx_conf(double *totplaq,quad_su3 *conf);

#endif
