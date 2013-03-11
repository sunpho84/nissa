#ifndef _WILSON_FORCE_H
#define _WILSON_FORCE_H

#include "../../new_types/new_types_definitions.h"

void Wilson_force_WORKER();
void Wilson_force(quad_su3 **F,quad_su3 **eo_conf,double beta);

#endif
