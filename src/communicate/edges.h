#ifndef _EDGES_H
#define _EDGES_H

#include "new_types/new_types_definitions.h"

void communicate_eo_edges(char **data,comm_t &comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
void communicate_eo_quad_su3_edges(quad_su3 **conf);
void communicate_lx_edges(char *data,comm_t &comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
void communicate_lx_quad_su3_edges(quad_su3 *conf);
void communicate_lx_su3_edges(su3 *u);
#endif
