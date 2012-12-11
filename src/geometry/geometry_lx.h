#ifndef _geometry_Ax_h
#define _geometry_Ax_h
#include <mpi.h>
#include "../new_types/new_types_definitions.h"
int bordlx_of_coord(int *x,int mu);
int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int mu);
int edgelx_of_coord(int *x,int mu,int nu);
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
void glb_coord_of_glblx(coords x,int gx);
void initialize_lx_bord_receivers_of_kind(MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base);
void initialize_lx_bord_senders_of_kind(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *base);
void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
void rank_coord_of_site_of_coord(coords rank_coord,coords glb_coord);
void set_lx_bord_senders_and_receivers(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base);
void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
void set_lx_geometry();
void unset_lx_geometry();
#endif
