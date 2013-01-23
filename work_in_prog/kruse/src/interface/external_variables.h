#ifndef _EXTERNAL_H
#define _EXTERNAL_H

#include "nissa_types.h"

extern int loc_size[4],glb_size[4];
extern int loc_vol,bord_vol;
extern int rank;
extern coords nrank_dir,rank_coord;
extern coords *glb_coord_of_loclx;
extern coords rank_neigh[2],rank_neighdw,rank_neighup;
extern _Complex double ka0,ka1,ka2,ka3;
//tobeadded
extern int *g_eo2lexic,*g_lexic2eo,*g_lexic2eosub;

#endif
