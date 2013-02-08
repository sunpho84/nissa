#ifndef _RENAMING_SCHEME_H
#define _RENAMING_SCHEME_H

#define  T (loc_size[0])
#define LX (loc_size[1])
#define LY (loc_size[2])
#define LZ (loc_size[3])
#define VOLUME loc_vol
#define RAND   bord_vol
#define VOLUMEPLUSRAND (loc_vol+bord_vol)
#define T_global (glb_size[0])
#define N_PROC_X (nrank_dir[1])
#define N_PROC_Y (nrank_dir[2])
#define N_PROC_Z (nrank_dir[3])
#define g_proc_coords rank_coord
#define g_coord glb_coord_of_loclx
#define g_lexic2eosub loceo_of_loclx
#define g_proc_id rank
#define g_nb_t_dn (rank_neighdw[0])
#define g_nb_x_dn (rank_neighdw[1])
#define g_nb_y_dn (rank_neighdw[2])
#define g_nb_z_dn (rank_neighdw[3])
#define g_nb_t_up (rank_neighup[0])
#define g_nb_x_up (rank_neighup[1])
#define g_nb_y_up (rank_neighup[2])
#define g_nb_z_up (rank_neighup[3])
#define g_cart_id cart_rank
#define g_cart_grid cart_comm

//types
#define su3 tmlQCD_su3

//routines
#define Index full_lx_of_coords_list

#endif
