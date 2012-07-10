#ifndef _BACKFIELD_H
#define _BACKFIELD_H

#include "../new_types/new_types_definitions.h"

void add_backfield_to_conf(quad_su3 **conf,quad_u1 **u1);
void init_backfield_to_id(quad_u1 **S);
void rem_backfield_from_conf(quad_su3 **conf,quad_u1 **u1);

#endif
