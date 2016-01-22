#ifndef _GEOMETRY_LX_HPP
#define _GEOMETRY_LX_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#include <stdint.h>

#ifndef EXTERN_GEOMETRY_LX
 #define EXTERN_GEOMETRY_LX extern
 #define ONLY_INSTANTIATION
#endif

namespace nissa
{
  typedef int coords[NDIM];
  typedef double momentum_t[NDIM];
  
  //nomenclature:
  //-glb is relative to the global grid
  //-loc to the local one
  EXTERN_GEOMETRY_LX coords glb_size,loc_size;
  EXTERN_GEOMETRY_LX int64_t glb_vol,glb_spat_vol,glb_volh;
  EXTERN_GEOMETRY_LX int64_t loc_vol,loc_spat_vol,loc_volh;
  EXTERN_GEOMETRY_LX int64_t bulk_vol,non_bw_surf_vol,non_fw_surf_vol;
  EXTERN_GEOMETRY_LX int64_t surf_vol,bw_surf_vol,fw_surf_vol;
  EXTERN_GEOMETRY_LX int64_t vsurf_vol,vsurf_volh;
  EXTERN_GEOMETRY_LX int64_t vdir_bord_vol,vdir_bord_volh;
  EXTERN_GEOMETRY_LX double glb_vol2,loc_vol2;
  //-lx is lexicografic
  EXTERN_GEOMETRY_LX coords *glb_coord_of_loclx;
  EXTERN_GEOMETRY_LX coords *loc_coord_of_loclx;
  EXTERN_GEOMETRY_LX int *glblx_of_loclx;
  EXTERN_GEOMETRY_LX int *glblx_of_bordlx;
  EXTERN_GEOMETRY_LX int *loclx_of_bordlx;
  EXTERN_GEOMETRY_LX int *surflx_of_bordlx;
  EXTERN_GEOMETRY_LX int *glblx_of_edgelx;
  EXTERN_GEOMETRY_LX int *loclx_of_bulklx;
  EXTERN_GEOMETRY_LX int *loclx_of_surflx;
  EXTERN_GEOMETRY_LX int *loclx_of_non_bw_surflx;
  EXTERN_GEOMETRY_LX int *loclx_of_non_fw_surflx;
  EXTERN_GEOMETRY_LX int *loclx_of_bw_surflx;
  EXTERN_GEOMETRY_LX int *loclx_of_fw_surflx;
  EXTERN_GEOMETRY_LX int lx_geom_inited;
  //box, division in 2^NDIM of the lattice
  EXTERN_GEOMETRY_LX coords box_coord[1<<NDIM];
  EXTERN_GEOMETRY_LX coords box_size[1<<NDIM];
  EXTERN_GEOMETRY_LX int nsite_per_box[1<<NDIM];
  //neighbours of local volume + borders
  EXTERN_GEOMETRY_LX coords *loclx_neighdw,*loclx_neighup;
  EXTERN_GEOMETRY_LX coords *loclx_neigh[2];
  //ranks
  EXTERN_GEOMETRY_LX coords fix_nranks;
  EXTERN_GEOMETRY_LX int rank,nranks,cart_rank;
  EXTERN_GEOMETRY_LX coords rank_coord;
  EXTERN_GEOMETRY_LX coords rank_neigh[2],rank_neighdw,rank_neighup;
  EXTERN_GEOMETRY_LX coords plan_rank,line_rank,line_coord_rank;
  EXTERN_GEOMETRY_LX coords nrank_dir;
  EXTERN_GEOMETRY_LX int grid_inited;
  EXTERN_GEOMETRY_LX int nparal_dir;
  EXTERN_GEOMETRY_LX coords paral_dir;
  //size of the border and edges
  EXTERN_GEOMETRY_LX int bord_vol,bord_volh;
  EXTERN_GEOMETRY_LX int edge_vol,edge_volh;
  //size along various dir
  EXTERN_GEOMETRY_LX int bord_dir_vol[NDIM],bord_offset[NDIM];
  EXTERN_GEOMETRY_LX int edge_dir_vol[NDIM*(NDIM+1)/2],edge_offset[NDIM*(NDIM+1)/2];
  EXTERN_GEOMETRY_LX int edge_numb[NDIM][NDIM];
  //mapping of ILDG data
  EXTERN_GEOMETRY_LX coords scidac_mapping;
  //perpendicular dir
  EXTERN_GEOMETRY_LX coords all_dirs;
#if NDIM >= 2
  EXTERN_GEOMETRY_LX int perp_dir[NDIM][NDIM-1];
#endif
#if NDIM >= 3
  EXTERN_GEOMETRY_LX int perp2_dir[NDIM][NDIM-1][NDIM-2];
#endif
#if NDIM >= 4
  EXTERN_GEOMETRY_LX int perp3_dir[NDIM][NDIM-1][NDIM-2][NDIM-3];
#endif
  EXTERN_GEOMETRY_LX int map_mu[4]
#ifndef ONLY_INSTANTIATION
  ={4,1,2,3}
#endif
    ;
  
  void get_stagphase_of_lx(coords ph,int ivol);
  int bordlx_of_coord(int *x,int mu);
  int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int mu);
  void coord_of_lx(coords x,int ilx,coords s);
  void coord_of_rank(coords c,int irank);
  int edgelx_of_coord(int *x,int mu,int nu);
  int full_lx_of_coords_list(const int t,const int x,const int y,const int z);
  int glblx_neighdw(int gx,int mu);
  int glblx_neighup(int gx,int mu);
  int glblx_of_comb(int b,int wb,int c,int wc);
  int glblx_of_coord(coords x);
  int glblx_of_coord_list(int x0,int x1,int x2,int x3);
  int glblx_of_diff(int b,int c);
  int glblx_of_summ(int b,int c);
  int glblx_opp(int b);
  int loclx_of_coord(coords x);
  int loclx_of_coord_list(int x0,int x1,int x2,int x3);
  int lx_of_coord(coords x,coords s);
  int rank_hosting_glblx(int gx);
  int rank_hosting_site_of_coord(coords x);
  int rank_of_coord(coords x);
  void get_loclx_and_rank_of_coord(int *ivol,int *rank,coords g);
  void get_loclx_and_rank_of_glblx(int *lx,int *rx,int gx);
  int get_glblx_of_rank_and_loclx(int irank,int loclx);
  void glb_coord_of_glblx(coords x,int gx);
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
  void rank_coord_of_site_of_coord(coords rank_coord,coords glb_coord);
  void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
  void set_lx_geometry();
  void unset_lx_geometry();
  void red_coords_of_hypercubic_red_point(coords h,int hyp_red);
  void lx_coords_of_hypercube_vertex(coords lx,int hyp_cube);
  int hypercubic_red_point_of_red_coords(coords h);
}

#undef EXTERN_GEOMETRY_LX

#endif
