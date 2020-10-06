#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "nissa.h"

#ifndef BGQ
 typedef complex bi_complex[2];
#else
 typedef vector4double bi_complex;
#endif

typedef bi_complex bi_color[3];
typedef bi_color bi_su3[3];
typedef bi_su3 bi_oct_su3[8];
typedef bi_color bi_spincolor[4];
typedef bi_color bi_halfspincolor[2];

#ifndef EXTERN
 #define EXTERN
#endif

#define REMAP_BARRIER 100328
#define HOPPING_MATRIX_APPLICATION_BARRIER 100329
 
EXTERN int *bgqlx_of_loclx,*loclx_of_bgqlx;
EXTERN bi_halfspincolor *bgq_hopping_matrix_output_binded;
EXTERN bi_halfspincolor **bgq_hopping_matrix_output_pointer;
EXTERN bi_halfspincolor *bgq_hopping_matrix_output_T_buffer;
EXTERN int bgqlx_t_vbord_vol,bgq_vsurf_vol;
