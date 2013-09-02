#ifndef _HOPPING_MATRIX_BGQ_EO_OR_OE_H
#define _HOPPING_MATRIX_BGQ_EO_OR_OE_H

#include "../new_types/new_types_definitions.h"

void apply_staggered_hopping_matrix_oe_or_eo_bgq_nocomm_nobarrier(bi_oct_su3 **conf,int istart,int iend,bi_color *in,int);
void bgq_staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill(bi_color *out,int);
void start_staggered_hopping_matrix_oe_or_eo_bgq_communications();
void finish_staggered_hopping_matrix_oe_or_eo_bgq_communications(int);

#endif
