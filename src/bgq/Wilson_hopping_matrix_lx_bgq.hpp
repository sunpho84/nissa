#ifndef _HOPPING_MATRIX_BGQ_LX_H
#define _HOPPING_MATRIX_BGQ_LX_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void apply_Wilson_hopping_matrix_lx_bgq_nocomm_nobarrier(bi_oct_su3 *conf,int istart,int iend,bi_spincolor *in);
  void bgq_Wilson_hopping_matrix_lx_vdir_VN_comm_and_buff_fill();
  void start_Wilson_hopping_matrix_lx_bgq_communications();
  void finish_Wilson_hopping_matrix_lx_bgq_communications();
}

#endif
