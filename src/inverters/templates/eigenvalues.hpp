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
    //scalar product of (in1,in2)
    void scalar_prod(complex out,complex *in1,complex *in2,complex *buffer,int loc_size)
    {
      GET_THREAD_ID();
      
      //local scalar product
      NISSA_PARALLEL_LOOP(i,0,loc_size)
	unsafe_complex_conj1_prod(buffer[i],in1[i],in2[i]);
      THREAD_BARRIER();
      
      //reduction
      complex_vector_glb_collapse(out,buffer,loc_size);
    }
    
    //Ai=Ai-Bi*c
    void complex_vector_subtssign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n)
    {
      GET_THREAD_ID();
      
      NISSA_PARALLEL_LOOP(i,0,n)
	complex_subt_the_prod(a[i],b[i],c);
      set_borders_invalid(a);
    }
    
    //orthogonalize v with respect to the nvec of V
    void modified_GS(complex *v,complex **V,complex *buffer,int nvec,int vec_size)
    {
      for(int i=0;i<nvec;i++)
	{
	  complex s;
	  scalar_prod(s,V[i],v,buffer,vec_size);
	  complex_vector_subtssign_complex_vector_prod_complex(v,V[i],s,vec_size);
	}
    }
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
    complex *buffer=nissa_malloc("buffer",mat_size, complex);
    
    //fill V and orthonormalize
    for(int i=0;i<wspace_max_size;i++)
      {
	filler(V[i]);
	modified_GS(V[i],V,buffer,i,mat_size);
	double rat;
	double_vector_normalize(&rat,(double*)(V[i]),(double*)(V[i]),1.0,2*mat_size);
      }
    
    //free workspace
    for(int i=0;i<wspace_max_size;i++) nissa_free(V[i]);
  }
}

#endif
