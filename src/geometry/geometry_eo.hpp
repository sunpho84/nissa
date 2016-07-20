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

#include "base/bench.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
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
  
  /////////////////////////////////////////////////// split and join ////////////////////////////////////////////
  
  //! do not compile for types of unmatched size
  template <class Tout,class Tin> void check_matching_size_eo()
  {static_assert(NBASE_EL(Tout)==NBASE_EL(Tin),"calling with different number of elements between types in and out");}
  
  //! paste the even and odd parts of a vector into a full lx vector
  template <class Tout,class Tin> THREADABLE_FUNCTION_2ARG(paste_eo_parts_into_lx_vector, Tout*,out_lx, Tin**,in_eo)
  {
    check_matching_size_eo<Tout,Tin>();
    GET_THREAD_ID();
    START_TIMING(remap_time,nremap);
    
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	for(int iel=0;iel<NBASE_EL(Tin);iel++)
	  ((FLATTENED_TYPE(Tout)*)out_lx)[loclx_of_loceo[par][ieo]][iel]=((FLATTENED_TYPE(Tin)**)in_eo)[par][ieo][iel];
    
    STOP_TIMING(remap_time);
    set_borders_invalid(out_lx);
  }
  THREADABLE_FUNCTION_END
  
  //! separate the even and odd part of a vector
  template <class Tout,class Tin> THREADABLE_FUNCTION_2ARG(split_lx_vector_into_eo_parts, Tout**,out_eo, Tin*,in_lx)
  {
    check_matching_size_eo<Tout,Tin>();
    GET_THREAD_ID();
    START_TIMING(remap_time,nremap);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int iel=0;iel<NBASE_EL(Tin);iel++)
	((FLATTENED_TYPE(Tout)**)out_eo)[loclx_parity[ivol]][loceo_of_loclx[ivol]][iel]=((FLATTENED_TYPE(Tin)*)in_lx)[ivol][iel];
    
    STOP_TIMING(remap_time);
    set_borders_invalid(out_eo);
  }
  THREADABLE_FUNCTION_END
  
  //! get the even or odd part of a vector
  template <class Tout,class Tin> THREADABLE_FUNCTION_3ARG(get_evn_or_odd_part_of_lx_vector, Tout*,out_eo, Tin*,in_lx, int,par)
  {
    check_matching_size_eo<Tout,Tin>();
    GET_THREAD_ID();
    START_TIMING(remap_time,nremap);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      if(loclx_parity[ivol]==par)
	for(int iel=0;iel<NBASE_EL(Tin);iel++)
	  ((FLATTENED_TYPE(Tout)*)out_eo)[loceo_of_loclx[ivol]][iel]=((FLATTENED_TYPE(Tin)*)in_lx)[ivol][iel];
    
    STOP_TIMING(remap_time);
    set_borders_invalid(out_eo);
  }
  THREADABLE_FUNCTION_END
}

#undef EXTERN_GEOMETRY_EO

#endif
