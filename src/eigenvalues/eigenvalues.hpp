#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_PARPACK
 #include "eigenvalues_autarchic.hpp"
#endif
#include "eigenvalues_arpack.hpp"

namespace nissa
{
#ifdef HAVE_ARPACK_PARPACK_H
#define NISSA_DEFAULT_USE_ARPACK 1
  extern int use_arpack;
#endif
  
  //use arpack or the autarchic implementation
  template <class Fmat,class Filler>
  void eigenvalues_of_hermatr_find(complex **eig_vec,double *eig_val,int neig,bool min_max,
					     const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
					     const double target_precision,const int niter_max,
					     const Filler &filler)
  {
    master_printf("Solving eigenproblem for %d %s eigenvalues,\n",neig,min_max?"max":"min");
    master_printf(" target precision: %lg\n",target_precision);
    
#ifdef USE_PARPACK
    if(use_arpack)
      {
	master_printf("Using arpack\n");
	eigenvalues_of_hermatr_find_arpack(eig_vec,eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,target_precision,niter_max,filler);
      }
    else
      {
#endif
	master_printf("Using autarchic implementation\n");
	eigenvalues_of_hermatr_find_autarchic(eig_vec,eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,target_precision,niter_max,filler);
	
	//close the scope
#ifdef USE_PARPACK
      }
#endif
  }
}

#endif
