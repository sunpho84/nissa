#ifndef _EDGES_H
#define _EDGES_H

#include "../new_types/new_types_definitions.h"

void communicate_eo_edges(char **data,MPI_Datatype *MPI_EO_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EO_BORDS_RECE,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
void communicate_eo_quad_su3_edges(quad_su3 **conf);
void communicate_lx_edges(char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
void communicate_lx_quad_su3_edges(quad_su3 *conf);
void communicate_lx_su3_edges(su3 *u);
#endif
