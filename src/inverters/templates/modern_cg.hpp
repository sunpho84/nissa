#ifndef _MODERN_CG_HPP
#define _MODERN_CG_HPP

#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //solve using conjugate gradient algorithm
  template <class Fmat> void cg_solve(complex *x,const Fmat& mat_impl,complex *b,int mat_size,int mat_size_to_allocate,double target_residue,int nmax_iters,int print_each)
  {
    crash("reimplement");
    // master_printf("Starting cg, target_residue: %lg, niter_max: %d\n",target_residue,nmax_iters);
    
    // complex *p=nissa_malloc("p",mat_size_to_allocate,complex);
    // complex *r=nissa_malloc("r",mat_size_to_allocate,complex);
    // complex *t=nissa_malloc("t",mat_size_to_allocate,complex);
    // complex *ap=nissa_malloc("ap",mat_size_to_allocate,complex);
    // vector_copy(r,b);
    // vector_copy(p,b);
    // vector_reset(x);
    
    // //compute norm of r and store it for computing relative residual
    // double rr=double_vector_glb_norm2(r,mat_size),rel_res;
    // double source_norm=sqrt(double_vector_glb_norm2(b,mat_size));
    // if(source_norm==0) crash("invalid norm");
    
    // //count iterations
    // int iter=1;
    // do
    //   {
    // 	//compute (p,Ap)
    // 	mat_impl(ap,p);
    // 	double pap;
    // 	double_vector_glb_scalar_prod(&pap,(double*)ap,(double*)p,mat_size*2);
	
    // 	//compute alpha, store rr as old one
    // 	double alpha=rr/pap;
    // 	double roro=rr;
	
    // 	//adjust new solution and residual
    // 	double_vector_summassign_double_vector_prod_double((double*)x,(double*)p,alpha,2*mat_size);
    // 	double_vector_summassign_double_vector_prod_double((double*)r,(double*)ap,-alpha,2*mat_size);
    // 	//compute new residual norm
    // 	rr=double_vector_glb_norm2(r,mat_size);
	
    // 	//adjust new krylov vector
    // 	double beta=rr/roro;
    // 	double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,beta,2*mat_size);
	
    // 	//compute relative residue
    // 	rel_res=sqrt((double)rr)/source_norm;
	
    // 	if(iter%print_each==0) master_printf("it: %d, res: %.16lg\n",iter,rel_res);
	
    // 	iter++;
    //   }
    // while(rel_res>target_residue and iter<nmax_iters);
    
    // nissa_free(ap);
    // nissa_free(r);
    // nissa_free(t);
    // nissa_free(p);
  }
}

#endif
