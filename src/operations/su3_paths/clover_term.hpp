#ifndef _CLOVER_TERM_HPP
#define _CLOVER_TERM_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  CUDA_HOST_AND_DEVICE void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,clover_term_t Cl,spincolor in);
  void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,clover_term_t *Cl,spincolor *in);
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,clover_term_t *Cl,colorspinspin *in);
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,clover_term_t *Cl,su3spinspin *in);
  void unsafe_apply_chromo_operator_to_spincolor_128(spincolor_128 *out,clover_term_t *Cl,spincolor_128 *in);
  void chromo_operator(clover_term_t *Cl,quad_su3 *conf);
  void chromo_operator(eo_ptr<clover_term_t> Cl,eo_ptr<quad_su3> conf);
  
  //include the factor cSW - note that we include the factor "-1/4" here
  inline double chromo_operator_factor(double cSW)
  {return -cSW/4;}
  inline void chromo_operator_adjust_cSW(clover_term_t *Cl,double cSW_new,double cSW_old)
  {/*master_printf("adjusting from: %lg to %lg\n",cSW_old,cSW_new);*/double_vector_prod_double((double*)Cl,(double*)Cl,chromo_operator_factor(cSW_new)/chromo_operator_factor(cSW_old),sizeof(clover_term_t)/sizeof(double)*loc_vol);}
  inline void chromo_operator_adjust_cSW(eo_ptr<clover_term_t> Cl,double cSW_new,double cSW_old)
  {/*master_printf("adjusting from: %lg to %lg\n",cSW_old,cSW_new);*/for(int eo=0;eo<2;eo++) double_vector_prod_double((double*)(Cl[eo]),(double*)(Cl[eo]),chromo_operator_factor(cSW_new)/chromo_operator_factor(cSW_old),sizeof(clover_term_t)/sizeof(double)*loc_volh);}
  inline void chromo_operator_include_cSW(clover_term_t *Cl,double cSW)
  {double_vector_prod_double((double*)Cl,(double*)Cl,-cSW/4,sizeof(clover_term_t)/sizeof(double)*loc_vol);}
  template <class T> void chromo_operator_include_cSW(T Cl,double cSW)
  {chromo_operator_adjust_cSW(Cl,cSW+1e-16,-4);}
  template <class T> void chromo_operator_remove_cSW(T Cl,double cSW)
  {chromo_operator_adjust_cSW(Cl,-4,cSW+1e-16);}
  
  template <class T1,class T2> void clover_term(T1 Cl,double cSW,T2 conf)
  {
    chromo_operator(Cl,conf);
    chromo_operator_include_cSW(Cl,cSW);
  }
  
  CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor(halfspincolor out,complex diag,clover_term_t Cl,halfspincolor in);
  CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor_128(halfspincolor_128 out,complex& diag,clover_term_t Cl,halfspincolor_128 in);
  CUDA_HOST_AND_DEVICE void unsafe_apply_point_chromo_operator_to_spincolor_128(spincolor_128 out,clover_term_t Cl,spincolor_128 in);
  CUDA_HOST_AND_DEVICE inline void apply_point_twisted_clover_term_to_halfspincolor(halfspincolor out,double mass,double kappa,clover_term_t Cl,halfspincolor in)
  {
    complex z={1/(2*kappa),mass};
    apply_point_diag_plus_clover_term_to_halfspincolor(out,z,Cl,in);
  }
  CUDA_HOST_AND_DEVICE inline void apply_point_twisted_clover_term_to_halfspincolor_128(halfspincolor_128 out,double mass,double kappa,clover_term_t Cl,halfspincolor_128 in)
  {
    complex z={1/(2*kappa),mass};
    apply_point_diag_plus_clover_term_to_halfspincolor_128(out,z,Cl,in);
  }
  
  CUDA_HOST_AND_DEVICE void invert_point_twisted_clover_term(inv_clover_term_t inv,double mass,double kappa,clover_term_t Cl);
  void invert_twisted_clover_term(inv_clover_term_t *inv,double mass,double kappa,clover_term_t *Cl);
}

#endif
