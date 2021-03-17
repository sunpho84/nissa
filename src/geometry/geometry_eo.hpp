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
  //Structure to hold an even/old field
  template <typename T>
  struct eo_ptr
  {
    //Type representing a pointer to type T
    using Tptr=T*;
    
    //Type representing a pair of pointers
    using Tptr2=Tptr[2];
    
    //Inner pointer pairs
    Tptr data[2];
    
    //Access to data[i]
    CUDA_HOST_AND_DEVICE Tptr& operator[](const int i)
    {
      static_assert(std::is_trivially_copyable<eo_ptr<T>>::value,"not trivially copyable");
      return data[i];
    }
    
    //Constant access to data[i]
    CUDA_HOST_AND_DEVICE const Tptr& operator[](const int i) const
    {
      return data[i];
    }
    
    //Create from a pair of pointers
    CUDA_HOST_AND_DEVICE eo_ptr(Tptr a,Tptr b) : data{a,b} {}
    
    //Create from an array of two pointers - to deprecate?
    // CUDA_HOST_AND_DEVICE eo_ptr(Tptr2 a)
    // {
    //   data[0]=a[0];
    //   data[1]=a[1];
    // }
    
    //Default creator
    CUDA_HOST_AND_DEVICE eo_ptr()
    {
    }
    
    //Check whether the two ptr are equals
    CUDA_HOST_AND_DEVICE bool operator==(const eo_ptr& oth) const
    {
      return oth[0]==data[0] and oth[1]==data[1];
    }
    
    //Check whether the two ptr are different
    CUDA_HOST_AND_DEVICE bool operator!=(const eo_ptr& oth) const
    {
      return not ((*this)==oth);
    }
  };
  
  //-eo is even-odd
  CUDA_MANAGED EXTERN_GEOMETRY_EO int *loclx_parity;
  CUDA_MANAGED EXTERN_GEOMETRY_EO int *loceo_of_loclx;
  CUDA_MANAGED EXTERN_GEOMETRY_EO eo_ptr<int> loclx_of_loceo;
  CUDA_MANAGED EXTERN_GEOMETRY_EO eo_ptr<int> surfeo_of_bordeo;
  CUDA_MANAGED EXTERN_GEOMETRY_EO eo_ptr<coords> loceo_neighup;
  CUDA_MANAGED EXTERN_GEOMETRY_EO eo_ptr<coords> loceo_neighdw;
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
