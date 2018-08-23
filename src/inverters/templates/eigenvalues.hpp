#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"

namespace nissa
{
  namespace internal_eigenvalues
  {
    void scalar_prod(complex out,complex *in1,complex *in2,complex *buffer,int loc_size);
    void complex_vector_subtssign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n);
    void modified_GS(complex *v,complex **V,complex *buffer,int nvec,int vec_size);
  }
  
  //find the neig eigenvalues closest to the target
  template <class Fmat,class Filler>
  void eigenvalues_herm_find(complex **eig_vec,double *eig_val,int neig,bool min_max,
			     const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
			     const double tol,const int niter_max,
			     const int wspace_min_size,const int wspace_max_size,
			     const Filler &filler)
  {
    using namespace internal_eigenvalues;
    
    //allocate workspace
    complex *V[wspace_max_size];
    for(int i=0;i<wspace_max_size;i++) V[i]=nissa_malloc("Vi",mat_size_to_allocate,complex);
    complex *buffer=nissa_malloc("buffer",mat_size,complex);
    complex *temp=nissa_malloc("temp",mat_size,complex);
    
    //fill V and orthonormalize
    double useless_rat;
    for(int i=0;i<wspace_max_size;i++)
      {
	filler(V[i]);
	modified_GS(V[i],V,buffer,i,mat_size);
	double_vector_normalize(&useless_rat,(double*)(V[i]),(double*)(V[i]),1.0,2*mat_size);
      }
    
    //generate interaction matrix
    int wspace_size=neig;
    complex M[wspace_size*wspace_size];
    for(int i=0;i<wspace_size;i++)
      {
	imp_mat(temp,V[i]);
	for(int j=0;j<=i;j++)
	  {
	    scalar_prod(M[i*wspace_size+j],V[j],temp,buffer,mat_size);
	    master_printf("(%lg,%lg)\t",M[i*wspace_size+j][RE],M[i*wspace_size+j][IM]);
	  }
	master_printf("\n");
      }
    
    //free workspace
    for(int i=0;i<wspace_max_size;i++) nissa_free(V[i]);
    nissa_free(buffer);
    nissa_free(temp);
  }
}

#endif
