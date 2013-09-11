#ifndef _GEOMETRY_EO_H
#define _GEOMETRY_EO_H

#include "new_types/new_types_definitions.h"

void addrem_stagphases_to_eo_conf(quad_su3 **eo_conf);
void filter_hypercube_origin_sites(color **vec);
void initialize_eo_bord_receivers_of_kind(MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base);
void initialize_eo_bord_senders_of_kind(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *base);
void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGES_RECE,MPI_Datatype *base);
void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base);
void set_eo_bord_senders_and_receivers(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base);
void set_eo_edge_senders_and_receivers(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *MPI_EO_EDGES_RECE,MPI_Datatype *base);
void set_eo_geometry();
void unset_eo_geometry();

#endif
