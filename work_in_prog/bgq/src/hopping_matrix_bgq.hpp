#ifndef _HOPPING_MATRIX_BGQ_H
#define _HOPPING_MATRIX_BGQ_H

void apply_Wilson_hppopping_matrix_bgq_binded_nocomm_nobarrier(bi_oct_su3 *conf,int istart,int iend,bi_spincolor *in);
void bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill(bi_hppalfspincolor *out);
void start_Wilson_hppopping_matrix_bgq_binded_communications();
void finish_Wilson_hppopping_matrix_bgq_binded_communications();

#endif
