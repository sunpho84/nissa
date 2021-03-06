#ifndef _GEOMETRY_LX_HPP
#define _GEOMETRY_LX_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#include <stdint.h>
#include <routines/math_routines.hpp>

#ifndef EXTERN_GEOMETRY_LX
 #define EXTERN_GEOMETRY_LX extern
 #define ONLY_INSTANTIATION
#endif

#define NISSA_LOC_VOL_LOOP(a) for(int a=0;a<locVol;a++)

namespace nissa
{
  template <typename T,
	    size_t N>
  struct my_array
  {
    T data[N];
    
    CUDA_HOST_AND_DEVICE inline T& operator[](const size_t i)
    {
      return data[i];
    }
    
    CUDA_HOST_AND_DEVICE const inline T& operator[](const size_t i) const 
    {
      return data[i];
    }
  };
  
  typedef my_array<bool,NDIM> which_dir_t;
  typedef my_array<int,NDIM> coords_t;
  typedef my_array<double,NDIM> momentum_t;
  
  //nomenclature:
  //-glb is relative to the global grid
  //-loc to the local one
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t glbSize,locSize;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t glbVol,glbSpatVol,glbVolh;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int64_t locVol,locSpatVol,locVolh;
  EXTERN_GEOMETRY_LX int64_t bulkVol,nonBwSurfVol,nonFwSurfVol;
  EXTERN_GEOMETRY_LX int64_t surfVol,bwSurfVol,fwSurfVol;
  EXTERN_GEOMETRY_LX double glb_vol2,loc_vol2;
  //-lx is lexicografic
  //box, division in 2^NDIM of the lattice
  EXTERN_GEOMETRY_LX coords_t box_coord[1<<NDIM];
  EXTERN_GEOMETRY_LX coords_t box_size[1<<NDIM];
  CUDA_MANAGED EXTERN_GEOMETRY_LX int nsite_per_box[1<<NDIM];
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t *glbCoordOfLoclx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t *locCoordOfLoclx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int *glblxOfLoclx;
  EXTERN_GEOMETRY_LX int *glblxOfBordlx;
  EXTERN_GEOMETRY_LX int *loclxOfBordlx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int *surflxOfBordlx;
  EXTERN_GEOMETRY_LX int *glblxOfEdgelx;
  EXTERN_GEOMETRY_LX int *loclxOfBulklx;
  EXTERN_GEOMETRY_LX int *loclxOfSurflx;
  EXTERN_GEOMETRY_LX int *loclxOfNonBwSurflx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int *loclxOfNonFwSurflx;
  EXTERN_GEOMETRY_LX int *loclxOfBwSurflx;
  CUDA_MANAGED EXTERN_GEOMETRY_LX int *loclxOfFwSurflx;
  EXTERN_GEOMETRY_LX int lxGeomInited;
  //neighbours of local volume + borders
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t *loclxNeighdw,*loclxNeighup;
  EXTERN_GEOMETRY_LX coords_t *loclx_neigh[2];
  //ranks
  EXTERN_GEOMETRY_LX coords_t fix_nranks;
  EXTERN_GEOMETRY_LX int rank,nranks,cart_rank;
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t rank_coord;
  EXTERN_GEOMETRY_LX coords_t rank_neigh[2],rank_neighdw,rank_neighup;
  EXTERN_GEOMETRY_LX coords_t plan_rank,line_rank,line_coord_rank;
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t nrank_dir;
  EXTERN_GEOMETRY_LX int grid_inited;
  EXTERN_GEOMETRY_LX int nparal_dir;
  EXTERN_GEOMETRY_LX coords_t paral_dir;
  //size of the border and edges
  EXTERN_GEOMETRY_LX int bord_vol,bord_volh;
  EXTERN_GEOMETRY_LX int edge_vol,edge_volh;
  //size along various dir
  EXTERN_GEOMETRY_LX int bord_dir_vol[NDIM],bord_offset[NDIM];
  EXTERN_GEOMETRY_LX int edge_dir_vol[NDIM*(NDIM+1)/2],edge_offset[NDIM*(NDIM+1)/2];
  CUDA_MANAGED EXTERN_GEOMETRY_LX int edge_numb[NDIM][NDIM];
  //mapping of ILDG data
  CUDA_MANAGED EXTERN_GEOMETRY_LX coords_t scidac_mapping;
  //perpendicular dir
  EXTERN_GEOMETRY_LX which_dir_t all_dirs;
  EXTERN_GEOMETRY_LX which_dir_t only_dir[NDIM];
  EXTERN_GEOMETRY_LX which_dir_t all_other_dirs[NDIM];
  EXTERN_GEOMETRY_LX which_dir_t all_other_spat_dirs[NDIM];
#if NDIM >= 2
  CUDA_MANAGED EXTERN_GEOMETRY_LX int perp_dir[NDIM][NDIM-1];
#endif
#if NDIM >= 3
  CUDA_MANAGED EXTERN_GEOMETRY_LX int perp2_dir[NDIM][NDIM-1][NDIM-2];
#endif
#if NDIM >= 4
  EXTERN_GEOMETRY_LX int perp3_dir[NDIM][NDIM-1][NDIM-2][NDIM-3];
#endif
  CUDA_MANAGED EXTERN_GEOMETRY_LX int igamma_of_mu[4]
#ifndef ONLY_INSTANTIATION
  ={4,1,2,3}
#endif
    ;
  
  CUDA_HOST_AND_DEVICE coords_t get_stagphase_of_lx(const int& ivol);
  CUDA_HOST_AND_DEVICE int get_stagphase_of_lx(const int& ivol,const int& mu);
  
  int bordlx_of_coord(const coords_t& x,const int& mu);
  int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int mu);
  coords_t coord_of_lx(int ilx,const coords_t s);
  coords_t coord_of_rank(const int& irank);
  
  inline coords_t coord_summ(const coords_t& a1,const coords_t& a2,const coords_t& l)
  {
    coords_t s;
    
    for(int mu=0;mu<NDIM;mu++) s[mu]=(a1[mu]+a2[mu])%l[mu];
    
    return s;
  }
  
  inline void coord_summassign(coords_t& s,const coords_t& a,const coords_t& l)
  {
    s=coord_summ(s,a,l);
  }
  
  int edgelx_of_coord(const coords_t& x,const int& mu,const int& nu);
  int full_lx_of_coords_list(const int t,const int x,const int y,const int z);
  int glblx_neighdw(const int& gx,const int& mu);
  int glblx_neighup(const int& gx,const int& mu);
  int glblx_of_comb(int b,int wb,int c,int wc);
  int glblx_of_coord(const coords_t& x);
  int glblx_of_coord_list(int x0,int x1,int x2,int x3);
  int glblx_of_diff(const int& b,const int& c);
  int glblx_of_summ(const int& b,const int& c);
  int glblx_opp(const int& b);
  CUDA_HOST_AND_DEVICE int loclx_of_coord(const coords_t& x);
  
  inline int loclx_of_coord_list(int x0,int x1,int x2,int x3)
  {
    const coords_t c={x0,x1,x2,x3};
    return loclx_of_coord(c);
  }
  
  CUDA_HOST_AND_DEVICE int lx_of_coord(const coords_t& x,const coords_t& s);
  int vol_of_lx(const coords_t& size);
  int rank_hosting_glblx(const int& gx);
  int rank_hosting_site_of_coord(const coords_t& x);
  int rank_of_coord(const coords_t& x);
  void get_loclx_and_rank_of_coord(int& ivol,int& rank,const coords_t& g);
  void get_loclx_and_rank_of_glblx(int& lx,int& rx,const int& gx);
  coords_t glb_coord_of_glblx(int gx);
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
  coords_t rank_coord_of_site_of_coord(const coords_t& glb_coord);
  void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void set_lx_geometry();
  void unset_lx_geometry();
  coords_t get_mirrorized_site_coords(const coords_t& c,const int& imir);
  coords_t red_coords_of_hypercubic_red_point(int hyp_red);
  coords_t lx_coords_of_hypercube_vertex(int hyp_cube);
  int hypercubic_red_point_of_red_coords(const coords_t& h);
  
  //get mirrorized coord
  inline int get_mirrorized_site_coord(const int& c,const int& mu,const bool& flip)
  {
    return (glbSize[mu]+(1-2*flip)*c)%glbSize[mu];
  }
  
  //get mirrorized coords according to a bit decomposition of imir
  inline coords_t get_mirrorized_site_coords(const coords_t& c,const int& imir)
  {
    coords_t cmir;
    
    for(int mu=0;mu<NDIM;mu++)
      cmir[mu]=get_mirrorized_site_coord(c[mu],mu,get_bit(imir,mu));
    
    return cmir;
  }
}

#undef EXTERN_GEOMETRY_LX

#endif
