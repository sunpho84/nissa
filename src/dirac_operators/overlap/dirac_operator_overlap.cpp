#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "dirac_operator_overlap_kernel_portable.hpp"
#include "dirac_operator_overlap_kernel2.hpp"
#include "linalgs/linalgs.hpp"

#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "communicate/borders.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "geometry/geometry_lx.hpp"
#include "inverters/overlap/cgm_invert_overlap_kernel2.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/remez/remez_algorithm.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

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
  // x^(-1/2) = a_0 + sum_k=1^n a_k (x + b_k)^(-1) where x \in [epsilon,1]
  //
  // Since the spectrum of H^2 is [lambda_min, lambda_max] we evaluate the approximation in [lambda_min/lambda_max,1]
  // and then we rescale the weights a_k by a factor 1/sqrt(lambda_max) to have an approximation valid in [lambda_min, lambda_max]
  //
  // Eventually, (H^2)^(-1/2) = a'_0 + sum_k=1^n a'_k (H^2 + b_k)^(-1) with a'_i = a_i / sqrt(lambda_max)

  THREADABLE_FUNCTION_5ARG(apply_overlap, spincolor*,out, quad_su3*,conf, double,M, double,minerr, spincolor*,in)
  {
    GET_THREAD_ID();
    
    communicate_lx_quad_su3_borders(conf);
    communicate_lx_spincolor_borders(in);
    
    complex lambda_min,lambda_max;
    rat_approx_t appr;
    int niter_max=1000000;
    double req_res=minerr;
    int mat_size=loc_vol*NCOL; //physical volume
    int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL; //volume to allocate
    
    complex *eigen_vector=nissa_malloc("eigen_vector",mat_size_to_allocate,complex); // here we save the eigen vector, although it is not used in this function
    
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    
    // Application of H^2 to a spincolor vector and then cast to a complex vector
    const auto imp_mat=[conf,temp,M](complex *out_lx,complex *in_lx)
      {
	apply_overlap_kernel2((spincolor*)out_lx,conf,M,(spincolor*)temp,0.0,(spincolor*)in_lx);
      };
    const auto filler=[](complex *out_lx){generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);};
    
    //lambda_min = min eigenvalue of H^2
    const double min_max_prec=1e-5;
    const int MIN=0,MAX=1;
    const int neig=1;
    eigenvalues_find(&eigen_vector,&lambda_min,neig,MIN,mat_size,mat_size_to_allocate,imp_mat,min_max_prec,niter_max,filler); //find lambda_min (min eigenvalue of H^2)
    
    //lambda_max = max eigenvalue of H^2
    eigenvalues_find(&eigen_vector,&lambda_max,neig,MAX,mat_size,mat_size_to_allocate,imp_mat,min_max_prec,niter_max,filler); //find lambda_max (max eigenvalue of H^2)
    nissa_free(eigen_vector);
    
    // Since H is hermitian, H^2 has only real eigenvalues, so we check that their imaginary part is 0
    double minimum=lambda_min[RE];
    master_printf("min eigenvalue (%lg,%lg)\n",lambda_min[RE],lambda_min[IM]);
    double maximum=lambda_max[RE];
    master_printf("max eigenvalue (%lg,%lg)\n",lambda_max[RE],lambda_max[IM]);
    
    int num=-1,den=2;
    const double remez_minmax_err=0.01;
    generate_approx(appr, minimum,maximum,num,den,minerr,remez_minmax_err); // we evaluate the rational approximation of 1/sqrt(x) in [epsilon,1]
    
    // sum the constant and all the shifts
    summ_src_and_all_inv_overlap_kernel2_cgm(temp,conf,M,&appr,niter_max,req_res,in);
    // z = a'_0 in + sum_k=1^n a'_k (H^2+b_k)^(-1) in ( in is the input vector and n = appr.degree() )
    
    // here we apply H to z, out = H z, thus now out = sign(H) in
    apply_overlap_kernel(out,conf,M,temp);
    
    // here we apply g5 to out and we add the input vector, thus now out = in + g_5 sign(H) in = ( Id + g5*sign(H) ) in = D_ov in
    // this is horrible but fast (cit. Sunpho)
    NISSA_PARALLEL_LOOP(X,0,loc_vol)
      for(int c=0;c<NCOL;c++)
	{
	  out[X][0][c][0]=+out[X][0][c][0]+in[X][0][c][0];
	  out[X][0][c][1]=+out[X][0][c][1]+in[X][0][c][1];
	  out[X][1][c][0]=+out[X][1][c][0]+in[X][1][c][0];
	  out[X][1][c][1]=+out[X][1][c][1]+in[X][1][c][1];
	  out[X][2][c][0]=-out[X][2][c][0]+in[X][2][c][0];
	  out[X][2][c][1]=-out[X][2][c][1]+in[X][2][c][1];
	  out[X][3][c][0]=-out[X][3][c][0]+in[X][3][c][0];
	  out[X][3][c][1]=-out[X][3][c][1]+in[X][3][c][1];
	}
    
    set_borders_invalid(out);
    
    nissa_free(temp);
  }
  THREADABLE_FUNCTION_END
}
