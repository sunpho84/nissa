#ifndef _HOPPING_MATRIX_BGQ_LX_HPP
#define _HOPPING_MATRIX_BGQ_LX_HPP

#include "new_types/su3.hpp"
#include "new_types/two_stage_computation.hpp"

#ifndef EXTERN_HOPPING_MATRIX_BGQ_LX
 #define EXTERN_HOPPING_MATRIX_BGQ_LX extern
#endif

namespace nissa
{
  EXTERN_HOPPING_MATRIX_BGQ_LX two_stage_computation_pos_t virlx_hopping_matrix_output_pos;
  EXTERN_HOPPING_MATRIX_BGQ_LX two_stage_computation_pos_t viroe_or_vireo_hopping_matrix_output_pos[2];
  
  void apply_Wilson_hopping_matrix_lx_bgq_nocomm(vir_oct_su3 *conf,int istart,int iend,vir_spincolor *in);
  void bgq_Wilson_hopping_matrix_lx_vdir_VN_comm_and_buff_fill();
  void start_Wilson_hopping_matrix_lx_bgq_communications();
  void finish_Wilson_hopping_matrix_lx_bgq_communications();
}

#undef EXTERN_HOPPING_MATRIX_BGQ_LX

#endif
