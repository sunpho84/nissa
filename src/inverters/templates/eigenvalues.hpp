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
    void complex_vector_summassign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n);
    inline void complex_vector_subtassign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n)
    {
      complex d={-c[RE],-c[IM]};
      complex_vector_summassign_complex_vector_prod_complex(a,b,d,n);
    }
    
    void modified_GS(complex *v,complex **V,complex *buffer,int nvec,int vec_size);
    void eigenvalues_of_hermatr_find_all_and_sort(double *lambda,complex *M,int size,double tau);
  }
  
  //find the neig eigenvalues closest to the target
  template <class Fmat,class Filler>
  void eigenvalues_of_hermatr_find(complex **eig_vec,double *eig_val,int neig,double tau,bool min_max,
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
    complex *residue=nissa_malloc("residue",mat_size,complex);
    complex *temp=nissa_malloc("temp",mat_size,complex);
    
    //fill V and orthonormalize
    double useless_rat;
    for(int i=0;i<wspace_max_size;i++)
      {
	filler(V[i]);
	modified_GS(V[i],V,buffer,i,mat_size);
	double_vector_normalize(&useless_rat,(double*)(V[i]),(double*)(V[i]),1.0,2*mat_size);
      }
    
    //main loop
    int neig_to_find=neig;
    do
      {
	//generate interaction matrix
	int wspace_size=neig;
	complex M[wspace_size*wspace_size];
	for(int i=0;i<wspace_size;i++)
	  {
	    imp_mat(temp,V[i]);
	    for(int j=0;j<=i;j++)
	      scalar_prod(M[j+wspace_size*i],V[j],temp,buffer,mat_size);
	  }
	
	//main loop
	int nconv=0;
	do
	  {
	    //find all eigenvalues of the reduced problem, sort them by distance with tau
	    double red_eig_val[wspace_size];
	    eigenvalues_of_hermatr_find_all_and_sort(red_eig_val,M,wspace_size,tau);
	    
	    //count of converged
	    int nconv_in_loop=0;
	    const int old_wspace_size=wspace_size;
	    
	    //combine
	    for(int i=0;i<neig_to_find;i++)
	      vector_reset(eig_vec[i+nconv]);
	    for(int j=0;j<old_wspace_size;j++)
	      for(int i=0;i<neig_to_find;i++)
		complex_vector_summassign_complex_vector_prod_complex(eig_vec[i+nconv],V[j],M[i+old_wspace_size*j],mat_size);
	    
	    for(int i=0;i<neig_to_find;i++)
	      {
		complex *e=eig_vec[i+nconv];
		imp_mat(residue,e);
		double_vector_summassign_double_vector_prod_double((double*)residue,(double*)e,-red_eig_val[i],mat_size*2);
		double residue_norm=sqrt(double_vector_glb_norm2(residue,mat_size));
		master_printf("i: %d, eig: %lg, res: %lg\n",i,red_eig_val[i],residue_norm);
	      }
	    //     /* Compute norm of the residual and update arrays convind/keepind*/
	    //     resnrm_old[act] = resnrm[act];
	    //     resnrm[act] = sqrt(square_norm((spinor*) r, N, 1));
	    //     if (resnrm[act] < tol)
	    //       {
	    // 	convind[conv] = act; 
	    // 	conv = conv + 1; 
	    //       }
	    // else{
	    //   keepind[keep] = act; 
	    //   keep = keep + 1; 
	    // }
	    
	    
	  }
	while(0);
      }
    while(0);
    
    //free workspace
    for(int i=0;i<wspace_max_size;i++) nissa_free(V[i]);
    nissa_free(buffer);
    nissa_free(temp);
  }
}

#endif
