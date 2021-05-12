#ifndef _GEOMETRY_EO_HPP
#define _GEOMETRY_EO_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#ifndef EXTERN_GEOMETRY_EO
# define EXTERN_GEOMETRY_EO extern
# define INIT_GEOMETRY_EO(ARGS...)
#else
# define INIT_GEOMETRY_EO(ARGS...) (ARGS)
#endif

#define NISSA_DEFAULT_USE_EO_GEOM 1

#define NISSA_LOC_VOLH_LOOP(a) for(LocEoSite a=0;a<locVolh;a++)

#include <geometry/geometry_lx.hpp>
#include <metaProgramming/nonConstMethod.hpp>
#include <new_types/su3.hpp>
#include <tensor/lookupTable.hpp>

namespace nissa
{
  DECLARE_COMPONENT(LocEoSite,int64_t,DYNAMIC);
  
  DECLARE_COMPONENT(BordEoSite,int64_t,DYNAMIC);
  
  DECLARE_COMPONENT(EdgeEoSite,int64_t,DYNAMIC);
  
  DECLARE_COMPONENT(Parity,int32_t,2);
  
  inline CUDA_DEVICE Parity EVN=0;
  
  inline CUDA_DEVICE Parity ODD=1;
  
#define FOR_BOTH_PARITIES(NAME,CORE...)		\
  FOR_ALL_COMPONENT_VALUES(Parity,NAME,CORE)
  
#define UNROLL_FOR_BOTH_PARITIES(NAME,CORE...)		\
  UNROLL_FOR_ALL_COMPONENT_VALUES(Parity,NAME,CORE)
  
  /// Half the local volume
  CUDA_MANAGED EXTERN_GEOMETRY_EO LocEoSite locVolh;
  
  /// Half the border volume
  EXTERN_GEOMETRY_EO BordEoSite bordVolh;
  
  /// "Half" the edge volume
  EXTERN_GEOMETRY_EO EdgeEoSite edgeVolh;
  
  /// Size of half the local volume exteneded with border
  EXTERN_GEOMETRY_EO LocEoSite locVolhWithBord;
  
  /// Size of half the local volume exteneded with border and edge
  EXTERN_GEOMETRY_EO LocEoSite locVolhWithBordAndEdge;
  
  /// Structure to hold an even/old field
  template <typename T>
  struct eo_ptr
  {
    /// Type representing a pointer to type T
    using Tptr=T*;
    
    /// Inner pointer pairs
    Tptr data[2];
    
    /// Constant access to data[i]
    CUDA_HOST_DEVICE
    const Tptr& operator[](const Parity& i) const
    {
      return data[i()];
    }
    
    /// Non-Constant access to data[i]
    CUDA_HOST_DEVICE
    Tptr& operator[](const Parity& i)
    {
      return data[i()];
    }
    
    /// Create from a pair of pointers
    CUDA_HOST_DEVICE eo_ptr(Tptr a,Tptr b) :
      data{a,b}
    {
    }
    
    /// Default creator
    CUDA_HOST_DEVICE eo_ptr()
    {
    }
    
    /// Check whether the two ptr are equals
    CUDA_HOST_DEVICE bool operator==(const eo_ptr& oth) const
    {
      return oth[0]==data[0] and oth[1]==data[1];
    }
    
    /// Check whether the two ptr are different
    CUDA_HOST_DEVICE bool operator!=(const eo_ptr& oth) const
    {
      return not ((*this)==oth);
    }
  };
  
  /// Structure to hold an even/old field
  template <typename Tc,
	    typename Fund>
  struct EoTensor
  {
    /// Type representing a pointer to type T
    using Tptr=Tensor<Tc,Fund>;
    
    /// Inner pointer pairs
    Tptr data[2];
    
    /// Constant access to data[i]
    CUDA_HOST_DEVICE
    const Tptr& operator[](const int i) const
    {
      return data[i];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(operator[],CUDA_HOST_DEVICE);
    
    /// Create from a pair of pointers
    template <typename...ITc>
    CUDA_HOST_DEVICE EoTensor(const TensorCompFeat<ITc>&...tc)
    {
      for(int eo=0;eo<2;eo++)
	data[eo].allocate(tc.deFeat()...);
    }
    
    /// Default creator
    CUDA_HOST_DEVICE EoTensor()
    {
    }
  };
    
  CUDA_MANAGED EXTERN_GEOMETRY_EO LookupTable<OfComps<LocLxSite>,Parity> loclx_parity;
  CUDA_MANAGED EXTERN_GEOMETRY_EO LookupTable<OfComps<LocLxSite>,LocEoSite> loceo_of_loclx;
  CUDA_MANAGED EXTERN_GEOMETRY_EO LookupTable<OfComps<Parity,LocEoSite>,LocLxSite> loclx_of_loceo;
  CUDA_MANAGED EXTERN_GEOMETRY_EO LookupTable<OfComps<Parity,BordEoSite>,LocEoSite> surfeo_of_bordeo;
  CUDA_MANAGED EXTERN_GEOMETRY_EO LookupTable<OfComps<Parity,LocEoSite,Dir>,LocEoSite> loceo_neighup;
  CUDA_MANAGED EXTERN_GEOMETRY_EO LookupTable<OfComps<Parity,LocEoSite,Dir>,LocEoSite> loceo_neighdw;
  
  EXTERN_GEOMETRY_EO bool eo_geom_inited;
  EXTERN_GEOMETRY_EO bool use_eo_geom;
  
  /// Return the border id given the local site beyond the local e/o volume
  INLINE_FUNCTION BordEoSite bordEoSiteOfExtendedLocEoSize(const LocEoSite& locEoSite)
  {
    return locEoSite()-locVolh();
  }
  
    /// Return the local site beyond the e/o local volume, corresponding to the passed border id
  INLINE_FUNCTION LocEoSite extenedLocEoSiteOfBordEoSite(const BordEoSite& bordEoSite)
  {
    return locVolh+bordEoSite();
  }
  
  void filter_hypercube_origin_sites(color **vec);
  Parity glblx_parity(const GlbLxSite& glx);
  Parity glb_coord_parity(const GlbCoords& c);
  void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGES_RECE,MPI_Datatype *base);
  void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base);
  void set_eo_edge_senders_and_receivers(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *MPI_EO_EDGES_RECE,MPI_Datatype *base);
  void set_eo_geometry();
  void unset_eo_geometry();
}

#undef EXTERN_GEOMETRY_EO

#endif
