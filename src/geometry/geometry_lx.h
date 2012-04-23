#pragma once

int bordlx_of_coord(int *x,int mu);
int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int mu);
int edgelx_of_coord(int *x,int mu,int nu);
int glblx_of_coord(coords x);
int glblx_of_coord_list(int x0,int x1,int x2,int x3);
int loclx_of_coord(coords x);
int loclx_of_coord_list(int x0,int x1,int x2,int x3);
int lx_of_coord(coords x,coords s);
int rank_hosting_site_of_coord(coords x);
int rank_of_coord(coords x);
void get_loclx_and_rank_of_coord(int *ivol,int *rank,coords g);
void initialize_lx_bord_receivers_of_kind(MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base);
void initialize_lx_bord_senders_of_kind(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *base);
void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
void rank_coord_of_site_of_coord(coords rank_coord,coords glb_coord);
void set_lx_bord_senders_and_receivers(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base);
void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
void set_lx_geometry();
void unset_lx_geometry();
