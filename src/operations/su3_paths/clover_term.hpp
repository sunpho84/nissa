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
  //apply the chromo operator to the passed spincolor
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_apply_point_chromo_operator_to_spincolor(A&& out,
						       const B& Cl,
						       const C& in)
  {
    unsafe_su3_prod_color(out[0],Cl[0],in[0]);
    su3_dag_summ_the_prod_color(out[0],Cl[1],in[1]);
    unsafe_su3_prod_color(out[1],Cl[1],in[0]);
    su3_subt_the_prod_color(out[1],Cl[0],in[1]);
    
    unsafe_su3_prod_color(out[2],Cl[2],in[2]);
    su3_dag_summ_the_prod_color(out[2],Cl[3],in[3]);
    unsafe_su3_prod_color(out[3],Cl[3],in[2]);
    su3_subt_the_prod_color(out[3],Cl[2],in[3]);
  }
  
  void unsafe_apply_chromo_operator_to_spincolor(LxField<spincolor>& out,
						 const LxField<clover_term_t>& Cl,
						 const LxField<spincolor>& in);
  
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
  
  template <typename O,
	    typename C,
	    typename I>
  CUDA_HOST_AND_DEVICE void apply_point_squared_twisted_clover_term_to_halfspincolor(O&& out,
										     const double mass,
										     const double kappa,
										     const C& Cl,
										     const I& in,
										     const int& offset)
  {
    spincolor temp;
    apply_point_twisted_clover_term_to_halfspincolor(temp,+mass,kappa,Cl,in,offset);
    apply_point_twisted_clover_term_to_halfspincolor(out ,-mass,kappa,Cl,temp,offset);
  }
  
  // CUDA_HOST_AND_DEVICE inline void apply_point_twisted_clover_term_to_halfspincolor_128(halfspincolor_128 out,double mass,double kappa,clover_term_t Cl,halfspincolor_128 in)
  // {
  //   complex z={1/(2*kappa),mass};
  //   apply_point_diag_plus_clover_term_to_halfspincolor_128(out,z,Cl,in);
  // }
  
  /// Form the inverse of the clover term
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void invert_point_twisted_clover_term(A&& inv,
					const double& mass,
					const double& kappa,
					const B& Cl)
  {
    for(int x_high_low=0;x_high_low<2;x_high_low++)
      for(int x_id=0;x_id<NDIRAC/2;x_id++)
	for(int x_ic=0;x_ic<NCOL;x_ic++)
	  {
	    //prepare the point source
	    halfspincolor b;
	    halfspincolor_put_to_zero(b);
	    b[x_id][x_ic][RE]=1;
	    double ori_rr=halfspincolor_norm2(b);
	    
	    //reset the solution
	    halfspincolor x;
	    halfspincolor_put_to_zero(x);
	    
	    //prepare r and p
	    halfspincolor r,p;
	    halfspincolor_copy(r,b);
	    halfspincolor_copy(p,b);
	    
	    //norm of r is 1
	    double rr=1;
	    
	    //count iterations
	    int iter=1;
	    [[ maybe_unused ]]
	    const int niter_max=200,niter_for_verbosity=20;
	    const double target_res=1e-32;
	    double res;
	    do
	      {
		//compute (p,Ap)
		halfspincolor ap;
		apply_point_squared_twisted_clover_term_to_halfspincolor(ap,mass,kappa,Cl,p,2*x_high_low);
		double pap=halfspincolor_scal_prod(p,ap);
		
		//compute alpha, store rr as old one
		double alpha=rr/pap;
		double roro=rr;
		
		//adjust new solution and residual,
		//compute new residual norm
		halfspincolor_summ_the_prod_double(x, p,+alpha);
		halfspincolor_summ_the_prod_double(r,ap,-alpha);
		rr=halfspincolor_norm2(r);
		
		//adjust new krylov vector
		double beta=rr/roro;
		halfspincolor_summ_the_prod_double(p,r,p,beta);
		
		//compute resiude
		res=rr/ori_rr;
		
		//write residual
		if(iter>=niter_for_verbosity)
#ifndef COMPILING_FOR_DEVICE
		  master_printf("iter %d rel residue: %lg\n",iter,res)
#endif
		    ;
		iter++;
	      }
	    while(res>=target_res and iter<niter_max);
	    if(iter>=niter_max)
#ifndef COMPILING_FOR_DEVICE
	      crash("exceeded maximal number of iterations %d, arrived to %d with residue %lg, target %lg",niter_max,iter,res,target_res);
 #else
	    __trap();
#endif
	    //halfspincolor ap;
	    //apply_point_squared_twisted_clover_term_to_halfspincolor(ap,mass,kappa,Cl+2*x_high_low,x);
	    //halfspincolor_summ_the_prod_double(ap,b,-1);
	    //double trr=halfspincolor_norm2(ap)/ori_rr;
	    //master_printf("true residue: %lg vs %lg\n",trr,rr);
	    
	    //copy the solution after removing the hermitian
	    halfspincolor temp;
	    apply_point_twisted_clover_term_to_halfspincolor(temp,-mass,kappa,Cl,x,2*x_high_low);
	    for(int id=0;id<NDIRAC/2;id++)
	      for(int ic=0;ic<NCOL;ic++)
		for(int ri=0;ri<2;ri++)
		  inv[x_high_low][id][ic][x_id][x_ic][ri]=temp[id][ic][ri];
	  }
  }
  
  template <typename InvClF,
	    typename ClF>
   void invert_twisted_clover_term(FieldFeat<InvClF>& invCl,
				   const double& mass,
				   const double& kappa,
				   const FieldFeat<ClF>& Cl)
   {
     verbosity_lv2_master_printf("Computing inverse clover term for quark of mass %lg and kappa %lg\n",mass,kappa);
     
     NISSA_PARALLEL_LOOP(X,0,Cl->nSites())
       invert_point_twisted_clover_term((*invCl)[X],mass,kappa,(*Cl)[X]);
     NISSA_PARALLEL_LOOP_END;
     
     invCl->invalidateHalo();
   }
}

#endif

