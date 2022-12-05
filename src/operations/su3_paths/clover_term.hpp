#ifndef _CLOVER_TERM_HPP
#define _CLOVER_TERM_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  CUDA_HOST_AND_DEVICE void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,clover_term_t Cl,spincolor in);
  void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,clover_term_t *Cl,spincolor *in);
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,clover_term_t *Cl,colorspinspin *in);
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,clover_term_t *Cl,su3spinspin *in);
  void unsafe_apply_chromo_operator_to_spincolor_128(spincolor_128 *out,clover_term_t *Cl,spincolor_128 *in);
  
  void chromo_operator(LxField<clover_term_t>& Cl,
		       const LxField<quad_su3>& conf);
  
  void chromo_operator(EoField<clover_term_t>& Cl_eo,
		       const EoField<quad_su3>& conf_eo);
  
  //include the factor cSW - note that we include the factor "-1/4" here
  inline double chromo_operator_factor(const double& cSW)
  {
    return -cSW/4;
  }
  
  inline void chromo_operator_adjust_cSW(LxField<clover_term_t>& Cl,
					 const double& cSW_new,
					 const double& cSW_old)
  {
    /*master_printf("adjusting from: %lg to %lg\n",cSW_old,cSW_new);*/
    Cl*=chromo_operator_factor(cSW_new)/chromo_operator_factor(cSW_old);
  }
  
  inline void chromo_operator_adjust_cSW(EoField<clover_term_t>& Cl,
					 const double& cSW_new,
					 const double& cSW_old)
  {
    /*master_printf("adjusting from: %lg to %lg\n",cSW_old,cSW_new);*/
    for(int eo=0;eo<2;eo++)
      Cl[eo]*=chromo_operator_factor(cSW_new)/chromo_operator_factor(cSW_old);
  }
  
  // inline void chromo_operator_include_cSW(clover_term_t *Cl,double cSW)
  // {double_vector_prod_double((double*)Cl,(double*)Cl,-cSW/4,sizeof(clover_term_t)/sizeof(double)*loc_vol);}
  template <typename T>
  void chromo_operator_include_cSW(T&& Cl,
				   const double& cSW)
  {
    chromo_operator_adjust_cSW(Cl,cSW+1e-16,-4);
  }
  
  template <typename T>
  void chromo_operator_remove_cSW(T&& Cl,
				  const double& cSW)
  {
    chromo_operator_adjust_cSW(Cl,-4,cSW+1e-16);
  }
  
  template <typename T1,
	    typename T2>
  void clover_term(T1&& Cl,
		   const double& cSW,
		   const T2& conf)
  {
    crash("reimplement");
    // chromo_operator(Cl,conf);
    // chromo_operator_include_cSW(Cl,cSW);
  }
  
  CUDA_HOST_AND_DEVICE void fill_point_twisted_clover_term(halfspincolor_halfspincolor out,int x_high_low,clover_term_t C,double mass,double kappa);
  CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor(halfspincolor out,complex diag,clover_term_t Cl,halfspincolor in);
  // CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor_128(halfspincolor_128 out,complex& diag,clover_term_t Cl,halfspincolor_128 in);
  // CUDA_HOST_AND_DEVICE void unsafe_apply_point_chromo_operator_to_spincolor_128(spincolor_128 out,clover_term_t Cl,spincolor_128 in);
  
  template <typename O,
	    typename C,
	    typename I>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void apply_point_twisted_clover_term_to_halfspincolor(O&& out,
							const double mass,
							const double kappa,
							const C& Cl,
							const I& in,
							const int& offset)
  {
    const complex z={1/(2*kappa),mass};
    
    auto o0=out[offset+0],o1=out[offset+1];
    auto i0=in[offset+0],i1=in[offset+1];
    auto c0=Cl[offset+0],c1=Cl[offset+1];
    
    unsafe_color_prod_complex(o0,i0,z);
    su3_summ_the_prod_color(o0,c0,i0);
    su3_dag_summ_the_prod_color(o0,c1,i1);
    
    unsafe_color_prod_complex(o1,i1,z);
    su3_summ_the_prod_color(o1,c1,i0);
    su3_subt_the_prod_color(o1,c0,i1);
  }
  
  // CUDA_HOST_AND_DEVICE inline void apply_point_twisted_clover_term_to_halfspincolor_128(halfspincolor_128 out,double mass,double kappa,clover_term_t Cl,halfspincolor_128 in)
  // {
  //   complex z={1/(2*kappa),mass};
  //   apply_point_diag_plus_clover_term_to_halfspincolor_128(out,z,Cl,in);
  // }
  
  CUDA_HOST_AND_DEVICE void invert_point_twisted_clover_term(inv_clover_term_t inv,double mass,double kappa,clover_term_t Cl);
  
   void invert_twisted_clover_term(EvnField<inv_clover_term_t>& invCl,
				   const double& mass,
				   const double& kappa,
				   const EvnField<clover_term_t>& Cl);
}

#endif
