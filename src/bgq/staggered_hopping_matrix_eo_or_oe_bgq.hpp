#ifndef _STAGGERED_HOPPING_MATRIX_BGQ_EO_OR_OE_HPP
#define _STAGGERED_HOPPING_MATRIX_BGQ_EO_OR_OE_HPP

#include "new_types/su3.hpp"

namespace nissa
{
#define DEFINE_STAGGERED_HOPPING_MATRIX(PREC,PREC_TYPE,VIR_32_64_OCT_SU3_T,VIR_32_64_COLOR_T) \
  void NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(VIR_32_64_OCT_SU3_T **conf,int istart,int iend,VIR_32_64_COLOR_T *in,int); \
  void NAME3(bgq,PREC,staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill)(); \
  void NAME3(start,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)(); \
  void NAME3(finish,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)(int); \
  void NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,staggered_D_subtract_from_mass2_times_in_bgq)(VIR_32_64_COLOR_T *out,PREC_TYPE mass2,VIR_32_64_COLOR_T *in); \
  void NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,staggered_D_bgq)(VIR_32_64_COLOR_T *out); \
  void NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,staggered_D_bgq)(VIR_32_64_COLOR_T *out,double coeff);
  
  DEFINE_STAGGERED_HOPPING_MATRIX(double,double,vir_oct_su3,vir_color);
  DEFINE_STAGGERED_HOPPING_MATRIX(single,float,vir_single_oct_su3,vir_single_color);
}

#endif
