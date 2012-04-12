#pragma once

//vectors parameters
#define nissa_vect_string_length 20
#define nissa_vect_alignment 16

//vector tags name
#define BORDERS_ALLOCATED 1
#define BORDERS_VALID 2
#define EDGES_ALLOCATED 4
#define EDGES_VALID 8
#define BORDERS_COMMUNICATED_AT_LEAST_ONCE 16

//ODD/EVN
#define EVN 0
#define ODD 1

//random number generator table length
#define ran2_ntab 32                                                                                                         

//sqrt(2)
#define RAD2 1.414213562373095048801688724209l

//wrapper to the internal routines
#define nissa_malloc(a,b,c) internal_nissa_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define nissa_free(a) a=internal_nissa_free(a,__FILE__,__LINE__)
#define crash(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define decript_MPI_error(...) internal_decript_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define master_printf(...) master_fprintf(stdout,__VA_ARGS__)
#define verbosity_master_printf(lv,...) if(nissa_verbosity>=lv) master_printf(__VA_ARGS__)
#define verbosity_lv1_master_printf(...) verbosity_master_printf(1,__VA_ARGS__)
#define verbosity_lv2_master_printf(...) verbosity_master_printf(2,__VA_ARGS__)
#define verbosity_lv3_master_printf(...) verbosity_master_printf(3,__VA_ARGS__)

#define if_nissa_vect_not_initialized() if(main_nissa_arr!=((char*)&main_nissa_vect)+sizeof(nissa_vect))

//loops
#define nissa_loc_volh_loop(a) for(int a=0;a<loc_volh;a++)
#define nissa_loc_vol_loop(a) for(int a=0;a<loc_vol;a++)

