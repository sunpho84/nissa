#ifndef _MACROS_H
#define _MACROS_H

#include <stdio.h>

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

//argument contactenation and naming
#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)
#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)

//vectors parameters
#define NISSA_VECT_STRING_LENGTH 20
#if defined(BGQ) && !defined(BGQ_EMU)
 #define NISSA_VECT_ALIGNMENT 64
#else
 #define NISSA_VECT_ALIGNMENT 16
#endif

//vector tags name
#define BORDERS_ALLOCATED 1
#define BORDERS_VALID 2
#define EDGES_ALLOCATED 4
#define EDGES_VALID 8
#define BORDERS_COMMUNICATED_AT_LEAST_ONCE 16
#define DO_NOT_SET_FLAGS 1
#define SEND_BACKWARD_BORD 1
#define SEND_FORWARD_BORD 2

//ILDG files
#define ILDG_MAGIC_NO                   0x456789ab
#define ILDG_MB_MASK                    ((uint16_t)0x80)
#define ILDG_ME_MASK                    ((uint16_t)0x40)

//ODD/EVN
#define EVN 0
#define ODD 1

//real/imag
#define RE 0
#define IM 1

//hmc
#define SEA_THEORY 0

//double per color, spincolor and quad_su3
#define NREALS_PER_COLOR 6
#define NREALS_PER_SPINCOLOR 24
#define NREALS_PER_QUAD_SU3 72

//random number generator table length
#define RAN2_NTAB 32

//communications benchmark
#ifdef BENCH
 #define START_COMMUNICATIONS_TIMING() do{if(IS_MASTER_THREAD) tot_comm_time-=take_time();}while(0)
 #define STOP_COMMUNICATIONS_TIMING() do{if(IS_MASTER_THREAD) tot_comm_time+=take_time();}while(0)
 #define GET_THREAD_ID_FOR_COMMUNICATIONS_TIMINGS() GET_THREAD_ID()
#else
 #define START_COMMUNICATIONS_TIMING()
 #define STOP_COMMUNICATIONS_TIMING()
 #define GET_THREAD_ID_FOR_COMMUNICATIONS_TIMINGS()
#endif

//constants
#define NISSA_DEFAULT_VERBOSITY_LV 1
#define NISSA_DEFAULT_USE_EO_GEOM 1
#define NISSA_DEFAULT_USE_128_BIT_PRECISION 0
#define NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED 1
#define NISSA_DEFAULT_WARN_IF_NOT_COMMUNICATED 0
#define NISSA_DEFAULT_USE_ASYNC_COMMUNICATIONS 1
#define NISSA_DEFAULT_VNODE_PARAL_DIR 0

#ifdef USE_VNODES
 #define NVNODES 2
#endif

//Pi
#ifndef M_PI
 #define M_PI           3.14159265358979323846
#endif

//Omelyan lambda
#define OMELYAN_LAMBDA 0.1931833

//sqrt(2)
#define RAD2 1.414213562373095048801688724209l

//avoid spilling
#define REORDER_BARRIER() __asm volatile ("")

//wrapper to the internal routines
#define nissa_malloc(a,b,c) (c*)internal_nissa_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define nissa_free(a) internal_nissa_free((char**)&(a),__FILE__,__LINE__)
#define crash(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define crash_printing_error(code,...) internal_crash_printing_error(__LINE__,__FILE__,code,__VA_ARGS__)
#define decript_MPI_error(...) internal_decript_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define master_printf(...) master_fprintf(stdout,__VA_ARGS__)

//add verbosity macro
#if MAX_VERBOSITY_LV>=1
#define VERBOSITY_LV1 (verbosity_lv>=1)
#else
 #define VERBOSITY_LV1 0
#endif
#if MAX_VERBOSITY_LV>=2
 #define VERBOSITY_LV2 (verbosity_lv>=2)
#else
 #define VERBOSITY_LV2 0
#endif
#if MAX_VERBOSITY_LV>=3
 #define VERBOSITY_LV3 (verbosity_lv>=3)
#else
 #define VERBOSITY_LV3 0
#endif

//wrappers for verbosity_lv?
#define verbosity_lv1_master_printf(...) do{if(VERBOSITY_LV1) master_printf(__VA_ARGS__);}while(0)
#define verbosity_lv2_master_printf(...) do{if(VERBOSITY_LV2) master_printf(__VA_ARGS__);}while(0)
#define verbosity_lv3_master_printf(...) do{if(VERBOSITY_LV3) master_printf(__VA_ARGS__);}while(0)

#define IF_VECT_NOT_INITIALIZED() if(main_arr!=((char*)&main_vect)+sizeof(nissa_vect))

//loops
#define NISSA_LOC_VOLH_LOOP(a) for(int a=0;a<loc_volh;a++)
#define NISSA_LOC_VOL_LOOP(a) for(int a=0;a<loc_vol;a++)
    
#define CRASH_IF_NOT_ALIGNED(a,b) if((long long int)(void*)a%b!=0) crash("alignement problem");

#endif
