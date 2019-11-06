#ifndef _EIGENVALUES_ALL_HPP
#define _EIGENVALUES_ALL_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_EIGEN
 #include <eigen3/Eigen/Dense>
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //default value of the workspace size
#define DEFAULT_EIGPROB_WSPACE_SIZE 100
  
  //fill a matrix using an implicit function F, and use Eigen to diagonalize it
  template <class F>
  void print_all_eigenstuff(const F &f,int mat_size)
  {
#ifndef USE_EIGEN
    
    crash("Need Eigen");
    
#else
    
    using namespace Eigen;
    
    if(nranks>1 or thread_pool_locked==false)
      crash("Cannot work in parallel");
    
    ComplexEigenSolver<MatrixXcd> solver;
    
    //fill the matrix to be diagonalized
    MatrixXcd matr(mat_size,mat_size);
    complex *test=nissa_malloc("test",mat_size,complex);
    complex *out=nissa_malloc("out",mat_size,complex);
    for(int i=0;i<mat_size;i++)
      {
	vector_reset(test);
	test[i][RE]=1.0;
	
	f(out,test);
	
	for(int j=0;j<mat_size;j++)
	  matr(j,i)=std::complex<double>(out[j][RE],out[j][IM]);
      }
    
    //diagonalize
    solver.compute(matr);
    
    //print eigenvectors
    for(int ieig=0;ieig<mat_size;ieig++)
      {
	master_printf("---\n");
	master_printf("%d\n",ieig);
	std::complex<double> lambda=solver.eigenvalues()(ieig);
	master_printf("%.16lg %.16lg\n",lambda.real(),lambda.imag());
	
	for(int i=0;i<mat_size;i++)
	  {
	    complex &c=test[i];
	    c[RE]=solver.eigenvectors()(i,ieig).real();
	    c[IM]=solver.eigenvectors()(i,ieig).imag();
	  }
	
	f(out,test);
	
	complex e;
	complex_vector_glb_scalar_prod(e,(complex*)test,(complex*)out,mat_size);
	
	complex_vector_subtassign_complex_vector_prod_complex((complex*)out,(complex*)test,e,mat_size);
	master_printf(" (%.16lg,%.16lg), residue: %lg\n",ieig,e[RE],e[IM],sqrt(double_vector_glb_norm2(out,mat_size)));
      }
#endif
  }
}

#endif
