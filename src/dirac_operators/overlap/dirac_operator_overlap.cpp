#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "dirac_operator_overlap_kernel_portable.hpp"
#include "dirac_operator_overlap_kernel2.hpp"
#include "linalgs/linalgs.hpp"

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
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

  THREADABLE_FUNCTION_4ARG(apply_overlap, spincolor*,out, quad_su3*,conf, double,M, double,minerr, spincolor*,in)
  {

    if(!check_borders_valid(conf)) communicate_lx_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_lx_spincolor_borders(in);

		double minimum;
		complex lambda_min, lambda_max;
		spincolor z;
		rat_approx_t appr;
		double niter_max=1.0e+08, req_res=minerr;
		complex *eigen_vector[1];
		int mat_size=loc_vol*NCOL; // physical volume
		int mat_size_to_allocate=(loc_vol+board_vol)*NCOL; // volume to allocate
		
		eigen_vector[0]=nissa_malloc("eigen_vector", mat_size_to_allocate, complex); // here we save the eigen vector, although it is not used in this function

		// Application of H^2 to a spincolor vector and then cast to a complex vector
		const auto imp_mat=[conf,M](complex *out_lx, complex *in_lx)
								{
									apply_overlap_kernel2((spincolor*)out_lx, conf, M, (spincolor*)in_lx);
									set_borders_invalid(out_lx);
								};
		const auto filler=[](complex *out_lx) {}
	
		// lambda_min = min eigenvalue of H^2	
		eigenvalues_of_hermatr_find(eigen_vector,&lambda_min,1,0,mat_size,mat_size_to_allocate,imp_mat,1.0e-05,1.0e+06,filler); // find lambda_min (min eigenvalue of H^2)
		
		// lambda_max = max eigenvalue of H^2
		eigenvalues_of_hermatr_find(eigen_vector,&lambda_max,1,1,mat_size,mat_size_to_allocate,imp_mat,1.0e-05,1.0e+06,filler); // find lambda_max (max eigenvalue of H^2)
		
		nissa_free(eigen_vector[0]);

		// Since H is hermitian, H^2 has only real eigenvalues, so we check that their imaginary part is 0	
		if( (lambda_min[IM]!=0) || (lambda_max[IM]!=0) ) exit(1);

		minimum = lambda_min[RE] / lambda_max[RE]; // here we define minimum = epsilon = Re(lambda_min) / Re(lambda_max) since eigenvalues of H^2 are real
		generate_approx(&appr, minimum, 1.0, -1, 2, minerr, 0.01); // we evaluate the rational approximation of 1/sqrt(x) in [epsilon,1]
		
		// here we rescale the weights of the rational approximation by a factor 1/sqrt(lambda_max) in order to get an approximation valid in [lambda_min,lambda_max]
		for(int i=0; i<appr.degree(); i++) appr.weights[i]/sqrt(lambda_max);		
 
		// sum the constant and all the shifts
		summ_src_and_all_inv_overlap_kernel2_cgm(z, conf, M, appr, niter_max, req_res, in); // capire cosa sono niter_max e req_res
		// z = a'_0 in + sum_k=1^n a'_k (H^2+b_k)^(-1) in ( in is the input vector and n = appr.degree() )

		// here we apply H to z, out = H z, thus now out = sign(H) in
		appy_overlap_kernel(out, conf, M, z);

		// here we apply g5 to out and we add the input vector, thus now out = in + g_5 sign(H) in = ( Id + g5*sign(H) ) in = D_ov in
		// this is horrible but fast (cit. Sunpho)
		GET_THREAD_ID();
		NISSA_PARALLEL_LOOP(X,0,loc_vol)
		{
		for(int c=0;c<3;c++)
    	{
      out[X][0][c][0] =  out[X][0][c][0] + in[X][0][c][0];
      out[X][0][c][1] =  out[X][0][c][1] + in[X][0][c][1];
      out[X][1][c][0] =  out[X][1][c][0] + in[X][1][c][0];
      out[X][1][c][1] =  out[X][1][c][1] + in[X][1][c][1];
      out[X][2][c][0] = -out[X][2][c][0] + in[X][2][c][0];
      out[X][2][c][1] = -out[X][2][c][1] + in[X][2][c][1];
      out[X][3][c][0] = -out[X][3][c][0] + in[X][3][c][0];
      out[X][3][c][1] = -out[X][3][c][1] + in[X][3][c][1];
    	}
		}
	set_borders_invalid(out);	
  THREADABLE_FUNCTION_END
	}
}
