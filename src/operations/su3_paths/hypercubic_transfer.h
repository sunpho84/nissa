#ifndef _HYPERCUBIC_TRANSFER_H
#define _HYPERCUBIC_TRANSFER_H

#include "new_types/new_types_definitions.h"

void get_covariant_transport_to_hypercube_origin(su3 path,coords c_hyp_ori,const coords c_hyp_red,quad_su3 **conf);
void get_covariant_transport_to_hypercube_origin(su3 path,int ivol,int hyp_red,quad_su3 **conf);
void gauge_transfer_in_hypercube_from_origin(color **out,quad_su3 **conf,int hyp_red,color **in);

#endif
