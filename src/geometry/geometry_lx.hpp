#ifndef _GEOMETRY_LX_HPP
#define _GEOMETRY_LX_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#ifdef USE_MPI
# include <mpi.h>
#endif

#include <stdint.h>

#include <new_types/coords.hpp>
#include <routines/math_routines.hpp>
#include <routines/mpi_routines.hpp>
#include <tensor/lookupTable.hpp>
#include <tensor/tensor.hpp>

#ifndef EXTERN_GEOMETRY_LX
# define EXTERN_GEOMETRY_LX extern
# define INIT_GEOMETRY_LX(ARGS...)
# define ONLY_INSTANTIATION
#else
# define INIT_GEOMETRY_LX(ARGS...) ARGS
#endif

#define NISSA_LOC_VOL_LOOP(a) for(LocLxSite a=0;a<locVol;a++)

#define FOR_ALL_GLB_TIMES(t) for(GlbCoord t=0;t<glbTimeSize;t++)

namespace nissa
{
  DECLARE_COMPONENT(GlbLxSite,int64_t,DYNAMIC);
  
  /// Global coordinates of a site
  using GlbCoords=Coords<GlbLxSite>;
  
  /// Global coordinate
  using GlbCoord=GlbLxSite;
  
  DECLARE_COMPONENT(LocLxSite,int64_t,DYNAMIC);
  
  /// Local coordinate
  using LocCoord=LocLxSite;
  
  /// Local coordinates of a site
  using LocCoords=Coords<LocLxSite>;
  
  DECLARE_COMPONENT(GlbEoSite,int64_t,DYNAMIC);
  
  DECLARE_COMPONENT(BordLxSite,int64_t,DYNAMIC);
  
  DECLARE_COMPONENT(EdgeLxSite,int64_t,DYNAMIC);
  
  DECLARE_COMPONENT(BulkLxSite,int64_t,DYNAMIC);
  
  //nomenclature:
  //-glb is relative to the global grid
  //-loc to the local one
  
  /// Global lattice hcube sizes
  CUDA_MANAGED EXTERN_GEOMETRY_LX GlbCoords glbSize INIT_GEOMETRY_LX({});
  
  /// Local lattice hcube sizes
  CUDA_MANAGED EXTERN_GEOMETRY_LX LocCoords locSize INIT_GEOMETRY_LX({});
  
  /// Global size in time direction
  inline const GlbCoord& glbTimeSize=
    glbSize(tDir);
  
  /// Local size in time direction
  inline const LocCoord& locTimeSize=
    locSize(tDir);
  
  /// Global 4D volume
  CUDA_MANAGED EXTERN_GEOMETRY_LX GlbLxSite glbVol;
  
  /// Global spatial volume
  EXTERN_GEOMETRY_LX GlbLxSite glbSpatVol;
  
  /// Local 4D volume
  CUDA_MANAGED EXTERN_GEOMETRY_LX LocLxSite locVol;
  
  /// Local spatial volume
  CUDA_MANAGED EXTERN_GEOMETRY_LX LocLxSite locSpatVol;
  
  /// Half the global volume
  EXTERN_GEOMETRY_LX GlbEoSite glbVolh;
  
  /// Bulk local volume
  EXTERN_GEOMETRY_LX BulkLxSite bulkVol;
  
  /// Number of sites not on the backward surface
  EXTERN_GEOMETRY_LX LocLxSite nonBwSurfVol;
  
  /// Number of sites not on the forward surface
  EXTERN_GEOMETRY_LX LocLxSite nonFwSurfVol;
  
  /// Number of site on the surface (not in the bulk)
  EXTERN_GEOMETRY_LX LocLxSite surfVol;
  
  /// Number of sites on the backward surface
  EXTERN_GEOMETRY_LX LocLxSite bwSurfVol;
  
  /// Number of sites on the forward surface
  EXTERN_GEOMETRY_LX LocLxSite fwSurfVol;
  
  /// Global coordinates of local sites
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite,Dir>,GlbCoord> glbCoordOfLoclx INIT_GEOMETRY_LX({});
  
  /// Local coordinates of local sites
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite,Dir>,LocCoord> locCoordOfLoclx INIT_GEOMETRY_LX({});
  
  /// Global site given the local site
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite>,GlbLxSite> glblxOfLoclx INIT_GEOMETRY_LX({});
  
  /// Global site given the border site
  EXTERN_GEOMETRY_LX Tensor<OfComps<BordLxSite>,GlbLxSite> glblxOfBordlx;
  
  /// Local site given the borders site
  EXTERN_GEOMETRY_LX Tensor<OfComps<BordLxSite>,LocLxSite> loclxOfBordlx;
  
  /// Local site adjacent to a given border site
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<BordLxSite>,LocLxSite> loclxSiteAdjacentToBordLx INIT_GEOMETRY_LX({});
  
  /// Global site corresponding to edge sites
  EXTERN_GEOMETRY_LX Tensor<OfComps<EdgeLxSite>,GlbLxSite> glblxOfEdgelx;
  
  /// Local site corresponding to bulk site
  EXTERN_GEOMETRY_LX Tensor<OfComps<BulkLxSite>,LocLxSite> loclxOfBulklx;
  
  /// Local site given a Non-Backward site on the surface
  EXTERN_GEOMETRY_LX Tensor<OfComps<LocLxSite>,LocLxSite> loclxOfNonBwSurflx;
  
  /// Local site given a Non-Forwward site on the surface
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite>,LocLxSite> loclxOfNonFwSurflx INIT_GEOMETRY_LX({});
  
  /// Local site given a Backward site on the surface
  EXTERN_GEOMETRY_LX Tensor<OfComps<LocLxSite>,LocLxSite> loclxOfBwSurflx;
  
  /// Local site given a Forwward site on the surface
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite>,LocLxSite> loclxOfFwSurflx INIT_GEOMETRY_LX({});
  
  /// Return the local site beyond the local volume, corresponding to the passed border id
  INLINE_FUNCTION LocLxSite extenedLocLxSiteOfBordLxSite(const BordLxSite& bordLxSite)
  {
    return locVol+bordLxSite();
  }
  
  /// Return the border id given the local site beyond the local volume
  INLINE_FUNCTION BordLxSite bordLxSiteOfExtendedLocLxSize(const LocLxSite& locLxSite)
  {
    return locLxSite()-locVol();
  }
  
  /// Up neighobour of a global site in the mu direction
  GlbLxSite glblxNeighup(const GlbLxSite& gx,const Dir& mu);
  
  /// Down neighobour of a global site in the mu direction
  GlbLxSite glblxNeighdw(const GlbLxSite& gx,const Dir& mu);
  
  /// Neighbours in the backward direction
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite,Dir>,LocLxSite> loclxNeighdw INIT_GEOMETRY_LX({});
  
  /// Neighbours in the forwkward direction
  CUDA_MANAGED EXTERN_GEOMETRY_LX LookupTable<OfComps<LocLxSite,Dir>,LocLxSite> loclxNeighup INIT_GEOMETRY_LX({});
  
  INLINE_FUNCTION CUDA_HOST_DEVICE
  const LookupTable<OfComps<LocLxSite,Dir>,LocLxSite>& loclxNeigh(const Orientation& ori) //nasty
  {
    const LookupTable<OfComps<LocLxSite,Dir>,LocLxSite>* ref[2]={&loclxNeighdw,&loclxNeighup};
     
     return *ref[ori.nastyConvert()];
  }
  
  /// Keep track of whether the grid is initialized
  EXTERN_GEOMETRY_LX bool gridInited;
  
  /// Keep track of whether the lexicographic geometry is initialized
  EXTERN_GEOMETRY_LX bool lxGeomInited;
  
  /// Number of parallelized directions
  EXTERN_GEOMETRY_LX Dir nparal_dir;
  
  EXTERN_GEOMETRY_LX Coords<bool> paral_dir;
  
  /// Size of the border
  EXTERN_GEOMETRY_LX BordLxSite bordVol;
  
  /// Size of the local volume exteneded with border
  EXTERN_GEOMETRY_LX LocLxSite locVolWithBord;
  
  /// Size of the local volume exteneded with border and edge
  EXTERN_GEOMETRY_LX LocLxSite locVolWithBordAndEdge;
  
  /// Return the edge id given the local site beyond the local volume and border
  INLINE_FUNCTION EdgeLxSite edgeLxSiteOfExtendedLocLxSize(const LocLxSite& locLxSite)
  {
    return locLxSite()-locVolWithBord();
  }
  
  /// Size of the edges
  EXTERN_GEOMETRY_LX EdgeLxSite edge_vol;
  
  /// Border volume along various direction
  EXTERN_GEOMETRY_LX Coords<BordLxSite> bord_dir_vol;
  
  /// Offset inside the border volume where the specific direction starts
  EXTERN_GEOMETRY_LX Coords<BordLxSite> bord_offset;
  
  EXTERN_GEOMETRY_LX LocLxSite edge_dir_vol[NDIM*(NDIM+1)/2],edge_offset[NDIM*(NDIM+1)/2];
  CUDA_MANAGED EXTERN_GEOMETRY_LX int edge_numb[NDIM][NDIM];
  
  /// Mapping of ILDG direction w.r.t native
  CUDA_MANAGED EXTERN_GEOMETRY_LX Coords<Dir> scidac_mapping INIT_GEOMETRY_LX({});
  
  //perpendicular dir
  EXTERN_GEOMETRY_LX Coords<bool> all_dirs;
  EXTERN_GEOMETRY_LX Coords<bool> only_dir[NDIM];
  EXTERN_GEOMETRY_LX Coords<bool> all_other_dirs[NDIM];
  EXTERN_GEOMETRY_LX Coords<bool> all_other_spat_dirs[NDIM];
  
#if NDIM >= 2
  CUDA_MANAGED EXTERN_GEOMETRY_LX int perp_dir[NDIM][NDIM-1];
#endif
  
#if NDIM >= 3
  CUDA_MANAGED EXTERN_GEOMETRY_LX int perp2_dir[NDIM][NDIM-1][NDIM-2];
#endif
  
#if NDIM >= 4
  EXTERN_GEOMETRY_LX int perp3_dir[NDIM][NDIM-1][NDIM-2][NDIM-3];
#endif
  
  constexpr Gamma igamma_of_mu(const Dir& dir)
  {
    switch(dir())
      {
      case 0:
	return 4;
	break;
      case 1:
	return 1;
	break;
      case 2:
	return 2;
	break;
      case 3:
	return 3;
	break;
    }
    return 0;
  }
  
  /// Return the staggered phases for a given site
  CUDA_HOST_DEVICE void get_stagphase_of_lx(Coords<int>& ph,const LocLxSite& ivol);
  
  /// Return the staggered phases for a given site
  CUDA_HOST_DEVICE int get_stagphase_of_lx(const LocLxSite& ivol,const Dir& mu);
  
  /// Return the index of site of coord x in the 3d space obtained projecting away mu
  LocLxSite spatLxOfProjectedCoords(const LocCoords& x,const Dir& mu);
  
  /// Given a lx site ilx, determine its coordinates in the box of size s
  template <typename LxSite>
  void coord_of_lx(Coords<LxSite>& x,LxSite /*Don't make it const ref*/ ilx,const Coords<LxSite>& s)
  {
    for(Dir mu=NDIM-1;mu>=0;mu--)
      {
	x(mu)=ilx%s(mu);
	ilx/=s(mu);
      }
  }
  
  void coord_of_rank(RankCoords& c,const Rank& x); //nasty
  
  /// Return the index of site of coord x in the 2d space obtained projecting away mu and nu
  template <typename LxSite>
  LxSite lineLxOfDoublyProjectedCoords(const Coords<LxSite>& x,const Dir& mu,const Dir& nu)
  {
    LxSite ilx=0;
    
    FOR_ALL_DIRS(rho)
      if(rho!=mu and rho!=nu)
	ilx=ilx*locSize(rho)+x(rho);
    
    return ilx;
  }
  
  /// Return the index of site of coord x in a box of sides s
  template <typename LxSite>
  CUDA_HOST_DEVICE constexpr
  LxSite lx_of_coord(const Coords<LxSite>& x,const Coords<LxSite>& s)
  {
    LxSite ilx=0;
    
    FOR_ALL_DIRS(mu)
      ilx=ilx*s(mu)+x(mu);
    
    return ilx;
  }
  
  /// Returns the global index of a site of coords x
  INLINE_FUNCTION
  GlbLxSite glblx_of_coord(const GlbCoords& x)
  {
    return lx_of_coord<GlbLxSite>(x,glbSize);
  }
  
  /// Returns the global index of a site of a list of coords
  INLINE_FUNCTION
  GlbLxSite glblx_of_coord_list(const GlbCoord& x0,const GlbCoord& x1,const GlbCoord& x2,const GlbCoord& x3)
  {
    GlbCoords c;
    
    c(Dir(0))=x0;
    c(xDir)=x1;
    c(yDir)=x2;
    c(zDir)=x3;
    
    return glblx_of_coord(c);
  }
  
  CUDA_HOST_DEVICE LocLxSite loclx_of_coord(const LocCoords& x);
  
  inline LocLxSite loclx_of_coord_list(const LocCoord& x0,const LocCoord& x1,const LocCoord& x2,const LocCoord& x3)
  {
    LocCoords c;
    
    c(Dir(0))=x0;
    c(xDir)=x1;
    c(yDir)=x2;
    c(zDir)=x3;
    
    return loclx_of_coord(c);
  }
  
  /// Return the volume of a given box
  template <typename I>
  I vol_of_lx(const Coords<I>& size)
  {
    I vol=1;
    
    FOR_ALL_DIRS(mu)
      vol*=size(mu);
    
    return vol;
  }
  
  Rank rank_hosting_glblx(const GlbLxSite& gx);
  Rank rank_hosting_site_of_coord(const GlbCoords& x);
  Rank rank_of_coord(const RankCoords& x);
  
  void get_loclx_and_rank_of_coord(LocLxSite& ivol,Rank& rank,const GlbCoords& g);
  
  /// Return the global index of site addressed by rank and loclx
  GlbLxSite get_glblx_of_rank_and_loclx(const int irank,const LocLxSite& loclx);
  
  /// Given a global lx site, returns its coordinates
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void glb_coord_of_glblx(GlbCoords& x,GlbLxSite /*Don't make it const reference*/ glbLx)
  {
    for(Dir mu=NDIM-1;mu>=0;mu--)
      {
	x(mu)=(glbLx()%glbSize(mu)());
	glbLx/=glbSize(mu)();
      }
  }
  
  /// Given a global lx site, returns the local coordinates and the rank hosting it
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void get_loclx_and_rank_of_glblx(LocLxSite& ivol,Rank& rank,const GlbLxSite& g)
  {
    GlbCoords c;
    glb_coord_of_glblx(c,g);
    get_loclx_and_rank_of_coord(ivol,rank,c);
  }
  
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
  void rank_coord_of_site_of_coord(RankCoords& rank_coord,const GlbCoords& glb_coord);
  void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void set_lx_geometry();
  void unset_lx_geometry();
  
  /// Gets mirrorized coord
  inline GlbCoord get_mirrorized_site_coord(const GlbCoord& c,const Dir& mu,const bool flip)
  { //nasty rename making it clear that refers to glb
    return (glbSize(mu)+(1-2*flip)*c)%glbSize(mu);
  }
  
  //get mirrorized coords according to a bit decomposition of imir
  inline void get_mirrorized_site_coords(GlbCoords& cmir,const GlbCoords& c,const int imir)
  {
    FOR_ALL_DIRS(mu)
      cmir(mu)=get_mirrorized_site_coord(c(mu),mu,get_bit(imir,mu()));
  }
}

#undef EXTERN_GEOMETRY_LX
#undef ONLY_INSTANTIATION
#undef INIT_GEOMETRY_LX

#endif
