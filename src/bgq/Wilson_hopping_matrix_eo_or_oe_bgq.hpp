#ifndef _WILSON_HOPPING_MATRIX_BGQ_EO_OR_OE_H
#define _WILSON_HOPPING_MATRIX_BGQ_EO_OR_OE_H

#include "new_types/new_types_definitions.hpp"
#include "base/macros.hpp"

namespace nissa
{
#define DEFINE_WILSON_HOPPING_MATRIX(PREC,PREC_TYPE,BI_32_64_OCT_SU3_T,BI_32_64_SPINCOLOR_T) \
  void NAME3(apply,PREC,Wilson_hopping_matrix_oe_or_eo_bgq_nocomm)(BI_32_64_OCT_SU3_T **conf,int istart,int iend,BI_32_64_SPINCOLOR_T *in,int oe_or_eo); \
  void NAME3(bgq,PREC,Wilson_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill)(); \
  void NAME3(start,PREC,Wilson_hopping_matrix_oe_or_eo_bgq_communications)(); \
  void NAME3(finish,PREC,Wilson_hopping_matrix_oe_or_eo_bgq_communications)(int); \
  void NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,Wilson_2D_bgq)(BI_32_64_SPINCOLOR_T *out); \
  
  DEFINE_WILSON_HOPPING_MATRIX(double,double,bi_oct_su3,bi_spincolor);
  DEFINE_WILSON_HOPPING_MATRIX(single,float,bi_single_oct_su3,bi_single_spincolor);
}

#endif
