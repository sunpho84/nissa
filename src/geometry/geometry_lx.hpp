#ifndef _GEOMETRY_LX_HPP
#define _GEOMETRY_LX_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <stdint.h>

#ifdef USE_MPI
# include <mpi.h>
#endif

#include <metaprogramming/globalVariable.hpp>

#include <routines/math_routines.hpp>

#ifndef EXTERN_GEOMETRY_LX
# define EXTERN_GEOMETRY_LX extern
# define ONLY_INSTANTIATION
#endif

#define NISSA_LOC_VOL_LOOP(a) for(int a=0;a<locVol;a++)

#define NDIM 4

namespace nissa
{
  /// Simple array needed to be used also on device
  template <typename T,
	    size_t N>
  struct MyArray
  {
    /// Internal data
    T data[N];
    
    /// Access to the i-th element
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    T& operator[](const size_t& i)
    {
      return data[i];
    }
    
    /// Constant access to the i-th element
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    const T& operator[](const size_t& i) const
    {
      return data[i];
    }
  };
  
  /// Type used to specify which direction to use
  typedef MyArray<bool,NDIM> WhichDirs;
  
  /// Coordinates
  typedef MyArray<int,NDIM> Coords;
  
  /// Momentum components
  typedef MyArray<double,NDIM> Momentum;
  
  /// Type to specify the number of sites in a box
  typedef MyArray<int,1<<NDIM> nsite_per_box_t;
  
  PROVIDE_GLOBAL_VAR(int64_t,locVol)
  
  PROVIDE_GLOBAL_VAR(int64_t,locSpatVol)
  
  PROVIDE_GLOBAL_VAR(int64_t,locVolh)
  
  PROVIDE_GLOBAL_VAR(Coords,locSize)
  
  PROVIDE_GLOBAL_VAR(int64_t,glbVol)
  
  PROVIDE_GLOBAL_VAR(int64_t,glbSpatVol)
  
  PROVIDE_GLOBAL_VAR(int64_t,glbVolh)
  
  PROVIDE_GLOBAL_VAR(Coords,glbSize)
  
  inline int64_t bulkVol;
  
  inline int64_t nonBwSurfVol;
  
  inline int64_t nonFwSurfVol;
  
  inline int64_t surfVol;
  
  inline int64_t bwSurfVol;
  
  inline int64_t fwSurfVol;
  
  inline double glbVol2;
  
  inline double locVol2;
  
  inline Coords boxCoord[1<<NDIM];
  
  inline Coords boxSize[1<<NDIM];
  
  PROVIDE_GLOBAL_VAR(nsite_per_box_t,nsite_per_box);
  
  CUDA_MANAGED EXTERN_GEOMETRY_LX Coords *glbCoordOfLoclx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX Coords *locCoordOfLoclx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *glblxOfLoclx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *glblxOfBordlx;
  EXTERN_GEOMETRY_LX int64_t *loclxOfBordlx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *surflxOfBordlx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *surflxOfEdgelx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *glblxOfEdgelx;
  EXTERN_GEOMETRY_LX int64_t *loclxOfBulklx;
  EXTERN_GEOMETRY_LX int64_t *loclxOfSurflx;
  EXTERN_GEOMETRY_LX int64_t *loclxOfNonBwSurflx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *loclxOfNonFwSurflx;
  EXTERN_GEOMETRY_LX int64_t *loclxOfBwSurflx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t *loclxOfFwSurflx;
  
  inline bool lxGeomInited;
  
  CUDA_MANAGED EXTERN_GEOMETRY_LX Coords *loclxNeighdw,*loclxNeighup;
  CUDA_MANAGED EXTERN_GEOMETRY_LX Coords *loclx_neigh[2];
  
  inline Coords fix_nranks;
  
  inline int rank;
  
  inline int nranks;
  
  inline int cartRank;
  
  inline Coords rankCoord;
  
  inline Coords rankNeigh[2];
  
  inline Coords rankNeighdw;
  
  inline Coords rankNeighup;
  
  PROVIDE_GLOBAL_VAR(Coords,nRanksDir);
  
  inline int gridInited;
  
  inline int nParalDir;
  
  PROVIDE_GLOBAL_VAR(Coords,isDirParallel);
  
  //size of the border and edges
  PROVIDE_GLOBAL_VAR(int64_t,bordVol);
  
  PROVIDE_GLOBAL_VAR(int64_t,bordVolh);
  
  PROVIDE_GLOBAL_VAR(int64_t,edgeVol);
  
  PROVIDE_GLOBAL_VAR(int64_t,edgeVolh);
  
  constexpr int nEdges=
	      NDIM*(NDIM-1)/2;
  
  inline int64_t bordDirVol[NDIM];
  
  inline int64_t bordOffset[NDIM];
  
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t edge_dir_vol[nEdges],edge_offset[nEdges];
  
  CUDA_MANAGED EXTERN_GEOMETRY_LX int edge_dirs[nEdges][2];
  CUDA_MANAGED EXTERN_GEOMETRY_LX bool isEdgeParallel[nEdges];
  EXTERN_GEOMETRY_LX int rank_edge_neigh[2][2][nEdges];
  CUDA_MANAGED EXTERN_GEOMETRY_LX int edge_numb[NDIM][NDIM];
  
  //mapping of ILDG data
  constexpr Coords scidacMapping{0,3,2,1};
  
  constexpr WhichDirs allDirs{1,1,1,1};
  
  constexpr WhichDirs onlyDir[NDIM]
    {{1,0,0,0},
     {0,1,0,0},
     {0,0,1,0},
     {0,0,0,1}};
  
  constexpr WhichDirs allOtherDirs[NDIM]
    {{0,1,1,1},
     {1,0,1,1},
     {1,1,0,1},
     {1,1,1,0}};
  
  constexpr WhichDirs allOtherSpatDirs[NDIM]
    {{0,1,1,1},
     {0,0,1,1},
     {0,1,0,1},
     {0,1,1,0}};
  
  using PerpDirs=
    MyArray<MyArray<int,NDIM-1>,NDIM>;
  
  PROVIDE_GLOBAL_VAR(PerpDirs,perpDirs);
  
  constexpr MyArray<MyArray<MyArray<int,NDIM-2>,NDIM-1>,NDIM> perp2Dirs
    {{{{{2,3},{1,3},{1,2}}},
      {{{2,3},{0,3},{0,2}}},
      {{{1,3},{0,3},{0,1}}},
      {{{1,2},{0,2},{0,1}}}}};
  
  constexpr MyArray<MyArray<MyArray<int,NDIM-2>,NDIM-1>,NDIM> perp3Dirs
    {{{{{3,2},{3,1},{2,1}}},
      {{{3,2},{3,0},{2,0}}},
      {{{3,1},{3,0},{1,0}}},
      {{{2,1},{2,0},{1,0}}}}};
  
  constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  int iGammaOfMu(const int& mu)
  {
    return
      Coords{4,1,2,3}[mu];
  }
  
  CUDA_HOST_AND_DEVICE Coords getStagphaseOfLx(const int64_t& ivol);
  
  CUDA_HOST_AND_DEVICE int getStagphaseOfLx(const int64_t& ivol,
					    const int& mu);
  
  int64_t bordlxOfCoord(const Coords& x,
			const int& mu);
  
  int bordlxOfCoordList(const int& x0,
			const int& x1,
			const int& x2,
			const int& x3,
			const int& mu);
  
  Coords coordOfLx(const int64_t& ilx,
		   const Coords& s);
  
  Coords coordOfRank(const int& irank);
  
  inline Coords coordsSum(const Coords& a1,
			  const Coords& a2,
			  const Coords& l)
  {
    Coords s;
    
    for(int mu=0;mu<NDIM;mu++)
      s[mu]=(a1[mu]+a2[mu])%l[mu];
    
    return s;
  }
  
  inline void coordsSummassign(Coords& s,
			       const Coords& a,
			       const Coords& l)
  {
    s=coordsSum(s,a,l);
  }
  
  int64_t edgelxOfCoord(const Coords &x,
			  const int &mu,
			  const int &nu);
  
  int fullLxOfCoordsList(const int t,
			 const int x,
			 const int y,
			 const int z);
  
  int64_t glblxNeighdw(const int64_t& gx,
			const int& mu);
  
  int64_t glblxNeighup(const int64_t& gx,
			const int& mu);
  
  int64_t glblxOfComb(const int64_t& b,
			const int& wb,
			const int64_t& c,
			const int& wc);
  
  int64_t glblxOfCoord(const Coords& x);
  
  int64_t glblxOfCoordList(const int& x0,
			      const int& x1,
			      const int& x2,
			      const int& x3);
  
  int64_t glblxOfDiff(const int64_t& b,
			const int64_t& c);
  
  int64_t glblxOfSum(const int64_t& b,
			const int64_t& c);
  
  int64_t glblxOpp(const int64_t& b);
  
  CUDA_HOST_AND_DEVICE int64_t loclxOfCoord(const Coords& x);
  
  inline int64_t loclxOfCoordList(const int& x0,
				  const int& x1,
				  const int& x2,
				  const int& x3)
  {
    return loclxOfCoord(Coords{x0,x1,x2,x3});
  }
  
  CUDA_HOST_AND_DEVICE int64_t lxOfCoord(const Coords& x,
					 const Coords& s);
  
  int64_t volOfLx(const Coords& size);
  
  int rankHostingGlblx(const int64_t& gx);
  
  int rankHostingSiteOfCoords(const Coords& x);
  
  int rankOfCoords(const Coords& x);
  
  Coords glbCoordOfGlblx(const int64_t& gx);
  
  /// Return the local site and rank containing the global coordinates
  inline std::pair<int,int64_t> getLoclxAndRankOfCoords(const Coords& g)
  {
    Coords l,p;
    for(int mu=0;mu<NDIM;mu++)
      {
	p[mu]=g[mu]/locSize[mu];
	l[mu]=g[mu]-p[mu]*locSize[mu];
      }
    
    const int rank=rankOfCoords(p);
    const int64_t ivol=loclxOfCoord(l);
    
    return {rank,ivol};
  }
  
  /// Return the local site and rank containing the global site
  inline std::pair<int,int64_t> getLoclxAndRankOfGlblx(const int64_t& gx)
  {
    return getLoclxAndRankOfCoords(glbCoordOfGlblx(gx));
  }
  
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
  Coords rankCoordsOfSiteOfCoord(const Coords& glb_coord);
  void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void set_lx_geometry();
  void unset_lx_geometry();
  Coords get_mirrorized_site_coords(const Coords& c,const int& imir);
  
  //get mirrorized coord
  inline int get_mirrorized_site_coord(const int& c,
				       const int& mu,
				       const bool& flip)
  {
    return (glbSize[mu]+(1-2*flip)*c)%glbSize[mu];
  }
  
  //get mirrorized coords according to a bit decomposition of imir
  inline Coords get_mirrorized_site_coords(const Coords& c,
					   const int& imir)
  {
    Coords cmir;
    
    for(int mu=0;mu<NDIM;mu++)
      cmir[mu]=get_mirrorized_site_coord(c[mu],mu,get_bit(imir,mu));
    
    return cmir;
  }
}

#undef EXTERN_GEOMETRY_LX

#endif
