#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef HAVE_ARPACK_PARPACK_H
 #include "eigenvalues_autarchic.hpp"
#endif
#include "eigenvalues_arpack.hpp"

namespace nissa
{
  template <class Fmat,class Filler>
  void eigenvalues_of_hermatr_find(complex **eig_vec,double *eig_val,int neig,bool min_max,
					     const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
					     const double target_precision,const int niter_max,
					     const Filler &filler)
  {
#ifdef HAVE_ARPACK_PARPACK_H
    eigenvalues_of_hermatr_find_arpack(eig_vec,eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,target_precision,niter_max,filler);
#else
    eigenvalues_of_hermatr_find_autarchic(eig_vec,eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,target_precision,niter_max,filler);
#endif
  }
}

#endif
