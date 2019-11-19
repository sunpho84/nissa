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
  //
  template <typename T>
  struct eo_ptr
  {
    using Tptr=T*;
    using Tptr2=Tptr[2];
    
    Tptr data[2];
    
    CUDA_HOST_AND_DEVICE Tptr& operator[](const int i)
    {
      static_assert(std::is_trivially_copyable<eo_ptr<T>>::value,"not trivially copyable");
      return data[i];
    }
    
    CUDA_HOST_AND_DEVICE const Tptr& operator[](const int i) const
    {
      return data[i];
    }
    
    eo_ptr(Tptr a,Tptr b) : data{a,b} {}
    
    eo_ptr(Tptr2 a)
    {
      data[0]=a[0];
      data[1]=a[1];
    }
    
    eo_ptr()
    {
    }
  };
  
  eo_ptr<double> a;
  
  //-eo is even-odd
  CUDA_MANAGED EXTERN_GEOMETRY_EO int *loclx_parity;
  CUDA_MANAGED EXTERN_GEOMETRY_EO int *loceo_of_loclx;
  CUDA_MANAGED EXTERN_GEOMETRY_EO int *loclx_of_loceo[2];
  CUDA_MANAGED EXTERN_GEOMETRY_EO int *surfeo_of_bordeo[2];
  CUDA_MANAGED EXTERN_GEOMETRY_EO coords *loceo_neighup[2];
  CUDA_MANAGED EXTERN_GEOMETRY_EO coords *loceo_neighdw[2];
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
