#ifndef _GEOMETRY_EO_HPP
#define _GEOMETRY_EO_HPP

#ifndef EXTERN_GEOMETRY_EO
 #define EXTERN_GEOMETRY_EO extern
#endif

//ODD/EVN
#define EVN 0
#define ODD 1

#define NISSA_DEFAULT_USE_EO_GEOM 1

#define NISSA_LOC_VOLH_LOOP(a) for(int a=0;a<loc_volh;a++)

#include "geometry_lx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //-eo is even-odd
  EXTERN_GEOMETRY_EO int *loclx_parity;
  EXTERN_GEOMETRY_EO int *loceo_of_loclx;
  EXTERN_GEOMETRY_EO int *loclx_of_loceo[2];
  EXTERN_GEOMETRY_EO int *surfeo_of_bordeo[2];
  EXTERN_GEOMETRY_EO coords *loceo_neighup[2];
  EXTERN_GEOMETRY_EO coords *loceo_neighdw[2];
  EXTERN_GEOMETRY_EO int eo_geom_inited;
  EXTERN_GEOMETRY_EO int use_eo_geom;
  
  void filter_hypercube_origin_sites(color **vec);
  int glblx_parity(int glx);
  int glb_coord_parity(coords c);
  void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGES_RECE,MPI_Datatype *base);
  void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base);
  void set_eo_edge_senders_and_receivers(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *MPI_EO_EDGES_RECE,MPI_Datatype *base);
  void set_eo_geometry();
  void unset_eo_geometry();
}

#undef EXTERN_GEOMETRY_EO

#endif
