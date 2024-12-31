#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "eigenvalues_all.hpp"

#ifdef USE_PARPACK
 #include "eigenvalues_parpack.hpp"
#endif
#include "eigenvalues_autarchic.hpp"

namespace nissa
{
#ifdef USE_PARPACK
 #define NISSA_DEFAULT_USE_PARPACK 1
#endif

  //use arpack or the autarchic implementation
  template <class Fmat,class Filler>
  void eigenvalues_find(complex **eig_vec,complex *eig_val,int neig,bool min_max,
				   const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
				   const double target_precision,const int niter_max,
				   const Filler &filler,int wspace_size=DEFAULT_EIGPROB_WSPACE_SIZE)
  {
    MASTER_PRINTF("Solving eigenproblem for %d %s eigenvalues,\n",neig,min_max?"max":"min");
    MASTER_PRINTF(" target precision: %lg\n",target_precision);
    
#ifdef USE_PARPACK
    if(use_parpack)
      {
	MASTER_PRINTF("Using parpack\n");
	eigenvalues_find_parpack(eig_vec,eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,target_precision,niter_max,filler,wspace_size);
      }
    else
      {
#endif
	MASTER_PRINTF("Using autarchic implementation\n");
	CRASH("Fix the algorithm! Plus, it works only for hermitean matrix");
	eigenvalues_find_autarchic(eig_vec,eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,target_precision,niter_max,filler);
	
	//close the scope
#ifdef USE_PARPACK
      }
#endif
  }
}

#endif
