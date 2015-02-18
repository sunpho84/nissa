#ifndef _GEOMETRY_EO_HPP
#define _GEOMETRY_EO_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void addrem_stagphases_to_eo_conf(quad_su3 **eo_conf);
  void filter_hypercube_origin_sites(color **vec);
  int glblx_parity(int glx);
  int glb_coord_parity(coords c);
  void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGES_RECE,MPI_Datatype *base);
  void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base);
  void set_eo_edge_senders_and_receivers(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *MPI_EO_EDGES_RECE,MPI_Datatype *base);
  void set_eo_geometry();
  void unset_eo_geometry();
}

#endif
