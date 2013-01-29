#ifndef _EXTERNAL_VARIABLES_H
#define _EXTERNAL_VARIABLES_H

#include "nissa_types.h"

#include <mpi.h>

extern int loc_size[4],glb_size[4];
extern int loc_vol,bord_vol;
extern int rank,cart_rank;
extern int *loceo_of_loclx,*loclx_parity;
extern coords nrank_dir,rank_coord;
extern coords *glb_coord_of_loclx;
extern coords rank_neigh[2],rank_neighdw,rank_neighup;
extern MPI_Comm cart_comm;

#endif
