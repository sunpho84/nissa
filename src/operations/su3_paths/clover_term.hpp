#ifndef _CLOVER_TERM_HPP
#define _CLOVER_TERM_HPP

#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,clover_term_t Cl,spincolor in);
  void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,clover_term_t *Cl,spincolor *in);
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,clover_term_t *Cl,colorspinspin *in);
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,clover_term_t *Cl,su3spinspin *in);
  void unsafe_apply_chromo_operator_to_spincolor_128(spincolor_128 *out,clover_term_t *Cl,spincolor_128 *in);
  void chromo_operator(clover_term_t *Cl,quad_su3 *conf);
  
  //include the factor cSW - note that we include the factor "-1/4" here
  inline void chromo_operator_include_cSW(clover_term_t *Cl,double cSW)
  {double_vector_prod_double((double*)Cl,(double*)Cl,-cSW/4,sizeof(clover_term_t)/sizeof(double)*loc_vol);}
  inline void chromo_operator_remove_cSW(clover_term_t *Cl,double cSW)
  {double_vector_prod_double((double*)Cl,(double*)Cl,-4/cSW,sizeof(clover_term_t)/sizeof(double)*loc_vol);}
  
  inline void clover_term(clover_term_t *Cl,double cSW,quad_su3 *conf)
  {
    chromo_operator(Cl,conf);
    chromo_operator_include_cSW(Cl,cSW);
  }
  
  void apply_point_diag_plus_clover_term_to_halfspincolor(halfspincolor out,complex diag,clover_term_t Cl,halfspincolor in);
  void unsafe_apply_point_chromo_operator_to_spincolor_128(spincolor_128 out,clover_term_t Cl,spincolor_128 in);
  inline void apply_point_twisted_clover_term_to_halfspincolor(halfspincolor out,double mass,double kappa,clover_term_t Cl,halfspincolor in)
  {
    complex z={1/(2*kappa),mass};
    apply_point_diag_plus_clover_term_to_halfspincolor(out,z,Cl,in);
  }
  
  void invert_point_twisted_clover_term(inv_clover_term_t inv,double mass,double kappa,clover_term_t Cl);
  void invert_twisted_clover_term(inv_clover_term_t *inv,double mass,double kappa,clover_term_t *Cl);
}

#endif
