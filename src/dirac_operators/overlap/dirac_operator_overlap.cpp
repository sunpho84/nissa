#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "dirac_operator_overlap_kernel_portable.hpp"
#include "dirac_operator_overlap_kernel2.hpp"
#include "linalgs/linalgs.hpp"

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "communicate/borders.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "geometry/geometry_lx.hpp"
#include "inverters/overlap/cgm_invert_overlap_kernel2.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  // Written by C. Bonanno and M. Cardinali
  //
  // Apply the overlap operator D_ov to a vector: out = D_ov in
  //
  // The overlap operator is defined as D_ov = Id + g5 sign(H)
  //
  // The overlap kernel operator is defined as H = g5 * (D_wilson - M*Id). Note that H^\dagger = H thus H^\dagger H = H^2
  //
  // The sign function is defined as sign(H) = H * (H^2)^(-1/2)
  //
  // We use a rational approximation of (H^2)^(-1/2) to evaluate the sign function:
  // x^(-1/2) = a_0 + sum_k=1^n a_k (x + b_k)^(-1) where x \in [lambda_min,lambda_max]
  // where lm, lM are the lowest and highest eigenvalues of H^2.
  //
  // Eventually, (H^2)^(-1/2) = a_0 + sum_k=1^n a_k (H^2 + b_k)^(-1)
  //
  // We used two different routines because parpack cannot support nested calls, and we can recycle the approximation
  // libraries that prevent innested calls
  void generate_rat_approx_for_overlap(quad_su3* conf,rat_approx_t*  appr,double mass_overlap,double maxerr)
  {
    communicate_lx_quad_su3_borders(conf);
    
    complex *lambda=nissa_malloc("lambda",2,complex);
    const int MIN=0,MAX=1;
    complex &lambda_min=lambda[MIN];
    complex &lambda_max=lambda[MAX];
    
    int niter_max=1000000;
    int mat_size=loc_vol*NCOL*NDIRAC; //physical volume
    int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL*NDIRAC; //volume to allocate
    
    complex *eigen_vector=nissa_malloc("eigen_vector",mat_size_to_allocate,complex); // here we save the eigen vector, although it is not used in this function
    
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    
    // Application of H^2 to a spincolor vector and then cast to a complex vector
    const auto imp_mat=[conf,temp,mass_overlap](complex *out_lx,complex *in_lx)
      {
	apply_overlap_kernel2((spincolor*)out_lx,conf,mass_overlap,(spincolor*)temp,0.0,(spincolor*)in_lx);
      };
    const auto filler=[](complex *out_lx){generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);};
    
    //lambda_min = min eigenvalue of H^2
    const int neig=1;
    eigenvalues_find(&eigen_vector,&lambda_min,neig,MIN,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler); //find lambda_min (min eigenvalue of H^2)
    
    //lambda_max = max eigenvalue of H^2
    eigenvalues_find(&eigen_vector,&lambda_max,neig,MAX,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler); //find lambda_max (max eigenvalue of H^2)
    
    nissa_free(eigen_vector);
    
    // Since H is hermitian, H^2 has only real eigenvalues, so we check that their imaginary part is 0
    double minimum=lambda_min[RE];
    master_printf("min eigenvalue (%lg,%lg)\n",lambda_min[RE],lambda_min[IM]);
    double maximum=lambda_max[RE];
    master_printf("max eigenvalue (%lg,%lg)\n",lambda_max[RE],lambda_max[IM]);
    
    int num=-1,den=2;
    generate_approx_of_maxerr(*appr,minimum,maximum,maxerr,num,den); // we evaluate the rational approximation of 1/sqrt(x) in [epsilon,1]
    
    nissa_free(temp);
    nissa_free(lambda);
  }
  
  //Verify the approximation, by applying twice the sign function
  void verify_rat_approx_for_overlap(quad_su3 *conf_lx,rat_approx_t &appr,double mass_overlap,double residue)
  {
    //generates the source and gets its norm
    spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
    generate_undiluted_source(in,RND_GAUSS,-1);
    double nin=double_vector_glb_norm2(in,loc_vol);
    
    //temporary and output
    spincolor *tmp=nissa_malloc("tmp",loc_vol+bord_vol,spincolor);
    spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
    
    ///apply twice the sign operator
    int niter_max=10000000;
    summ_src_and_all_inv_overlap_kernel2_cgm(tmp,conf_lx,mass_overlap,&appr,niter_max,residue,in);
    apply_overlap_kernel(out,conf_lx,mass_overlap,tmp);
    summ_src_and_all_inv_overlap_kernel2_cgm(tmp,conf_lx,mass_overlap,&appr,niter_max,residue,out);
    apply_overlap_kernel(out,conf_lx,mass_overlap,tmp);
    
    //subtracts the results from the source and gets its norm
    double_vector_subtassign((double*)out,(double*)in,sizeof(spincolor)/sizeof(double)*loc_vol);
    double nout=double_vector_glb_norm2(out,loc_vol);
    
    master_printf("Norm of the source: %.16lg\n",sqrt(nin));
    master_printf("Norm of the difference: %.16lg\n",sqrt(nout));
    master_printf("Relative norm of the difference: %.16lg\n",sqrt(nout/nin));
    
    nissa_free(in);
    nissa_free(tmp);
    nissa_free(out);
  }
  
  void apply_overlap(spincolor* out,quad_su3* conf,rat_approx_t*  appr,double req_res,double mass_overlap,double mass,spincolor* in)
  {
    
    communicate_lx_quad_su3_borders(conf);
    communicate_lx_spincolor_borders(in);
    
    int niter_max=1000000;
    
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    
    // sum the constant and all the shifts
    summ_src_and_all_inv_overlap_kernel2_cgm(temp,conf,mass_overlap,appr,niter_max,req_res,in);
    // temp = a_0 in + sum_k=1^n a_k (H^2+b_k)^(-1) in ( in is the input vector and n = appr.degree() )
    // here we apply H to temp, out = H temp, thus now out = sign(H) in
    apply_overlap_kernel(out,conf,mass_overlap,temp);
    
    // here we apply g5 to out and we add the input vector, thus now out = in + g_5 sign(H) in = ( Id + g5*sign(H) ) in = D_ov in
    // this is horrible but fast (cit. Sunpho)
    NISSA_PARALLEL_LOOP(X,0,loc_vol)
      for(int c=0;c<NCOL;c++)
	{
	  out[X][0][c][0]=+out[X][0][c][0]+(1.0+mass)*in[X][0][c][0];
	  out[X][0][c][1]=+out[X][0][c][1]+(1.0+mass)*in[X][0][c][1];
	  out[X][1][c][0]=+out[X][1][c][0]+(1.0+mass)*in[X][1][c][0];
	  out[X][1][c][1]=+out[X][1][c][1]+(1.0+mass)*in[X][1][c][1];
	  out[X][2][c][0]=-out[X][2][c][0]+(1.0+mass)*in[X][2][c][0];
	  out[X][2][c][1]=-out[X][2][c][1]+(1.0+mass)*in[X][2][c][1];
	  out[X][3][c][0]=-out[X][3][c][0]+(1.0+mass)*in[X][3][c][0];
	  out[X][3][c][1]=-out[X][3][c][1]+(1.0+mass)*in[X][3][c][1];
	}
    NISSA_PARALLEL_LOOP_END;
    
    master_printf("Diagonal part of overlap operator: %lg\n",(1.0+mass));
    
    set_borders_invalid(out);
    
    nissa_free(temp);
  }
}
