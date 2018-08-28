#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"

#include "modern_cg.hpp"

// #ifdef HAVE_EIGEN_DENSE
 #include "eigen3/Eigen/Dense"
// #endif


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
    double iterated_classical_GS(complex *v,int vec_size,int nvec,complex **A,complex* buffer,const int max_cgs_it);
    void eigenvalues_of_hermatr_find_all_and_sort(complex *eig_vec,double *lambda,const complex *M,const int M_size,const int neig,const double tau);
    void combine_basis_to_restart(int nout,int nin,complex *coeffs,complex **vect,int vec_length);
    
    template <class Fmat>
    void all_eigenvalues_finder(const int mat_size,const Fmat &imp_mat)
    {
      using namespace Eigen;
      SelfAdjointEigenSolver<MatrixXcd> solver;
      
      complex *out=nissa_malloc("out",mat_size,complex);
      complex *in=nissa_malloc("in",mat_size,complex);
      
      //fill the matrix to be diagonalized
      MatrixXcd matr(mat_size,mat_size);
      for(int j=0;j<mat_size;j++)
	{
	  vector_reset(in);
	  in[j][RE]=1.0;
	  imp_mat(out,in);
	  for(int i=0;i<mat_size;i++)
	    matr(i,j)=std::complex<double>(out[i][RE],out[i][IM]);
	}
      
      std::cout<<matr.topLeftCorner(12,12)<<std::endl;
      
      //diagonalize
      solver.compute(matr);
      
      std::cout<<solver.eigenvalues()<<std::endl;
      
      nissa_free(out);
      nissa_free(in);
    }
  }
  
  //find the neig eigenvalues closest to the target
  template <class Fmat,class Filler>
  void eigenvalues_of_hermatr_find(complex **eig_vec,double *eig_val,int neig,double tau,bool min_max,
				   const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
				   const double target_precision,const int niter_max,const int linit_max,
				   const double toldecay,const double eps_tr,
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
    for(int i=0;i<wspace_max_size;i++)
      {
	//fill
	filler(V[i]);
	//orthogonalize
	modified_GS(V[i],V,buffer,i,mat_size);
	//normalize
	double useless_rat;
	double_vector_normalize(&useless_rat,(double*)(V[i]),(double*)(V[i]),1.0,2*mat_size);
      }
    
    //generate interaction matrix
    int wspace_size=wspace_min_size;
    complex M[wspace_max_size*wspace_max_size];
    for(int j=0;j<wspace_size;j++)
      {
	imp_mat(temp,V[j]);
	for(int i=j;i<wspace_size;i++)
	  {
	    scalar_prod(M[j+wspace_max_size*i],V[i],temp,buffer,mat_size);
	    complex_conj(M[i+wspace_max_size*j],M[j+wspace_max_size*i]);
	  }
      }
    
    //main loop
    int solvestep=1;
    int iter=0;
    int neig_conv=0;
    do
      {
	master_printf("Iteration %d, wspace size: %d [%d:%d]\n",iter,wspace_size,wspace_min_size,wspace_max_size);
	
	//reset residue norm
	double residue_norm=0.0;
	
	//find all eigenvalues of the reduced problem, sort them by distance with tau
	double red_eig_val[wspace_size];
	complex red_eig_vec[wspace_size*wspace_size];
	eigenvalues_of_hermatr_find_all_and_sort(red_eig_vec,red_eig_val,M,wspace_max_size,wspace_size,tau);
	
	//combine the vectors
	complex *e=eig_vec[neig_conv];
	vector_reset(e);
	for(int j=0;j<wspace_size;j++)
	  complex_vector_summassign_complex_vector_prod_complex(e,V[j],red_eig_vec[wspace_size*j],mat_size);
	
	//compute the residue
	double residue_norm_old=residue_norm;
	imp_mat(residue,e);
	double_vector_summassign_double_vector_prod_double((double*)residue,(double*)e,-red_eig_val[0],mat_size*2);
	residue_norm=sqrt(double_vector_glb_norm2(residue,mat_size));
	master_printf("eig: %lg, res: %lg\n",red_eig_val[0],residue_norm);
	
	//if converged
	if(residue_norm<target_precision)
	  {
	    //store eigenvalue
	    eig_val[neig_conv]=red_eig_val[0];
	    master_printf("Eigenvalue %d/%d, %lg converged!\n",neig_conv,neig,eig_val[neig_conv]);
	    
	    //shift the others back
	    for(int i=0;i<wspace_size-1;i++)
	      red_eig_val[i]=red_eig_val[i+1];
	    
	    //restart using the remaining eigenvectors as basis
	    combine_basis_to_restart(wspace_size-1,wspace_size,red_eig_vec+1,V,mat_size);
	    
	    //update workspace size
	    wspace_size--;
	    
	    //make M diagonal
	    memset(M,0,sizeof(complex)*wspace_max_size*wspace_max_size);
	    for(int i=0;i<wspace_size;i++)
	      complex_put_to_real(M[i+i*wspace_max_size],red_eig_val[i]);
	    
	    //set tau closer to the minimal eig_val
	    if(min_max==0)
	      tau=std::max(eig_val[neig_conv],tau);
	    else
	      tau=std::min(eig_val[neig_conv],tau);
	    
	    //reset the stopping criterion
	    solvestep=1;
	    
	    neig_conv++;
	  }
	
	//reset if exceeded the workspace size
	if(wspace_size==wspace_max_size)
	  {
	    master_printf("Resetting\n");
	    
	    // master_printf("//////////////////////////// bef matr //////////////////////////\n");
	    // for(int i=0;i<wspace_size;i++)
	    // 	{
	    // 	  for(int j=0;j<=i;j++)
	    // 	      master_printf("(%+4.4g,%+4.4g)\t",M[j+wspace_max_size*i][RE],M[j+wspace_max_size*i][IM]);
	    // 	  master_printf("\n");
	    // 	}
	    
	    // master_printf("//////////////////////////// diag matr //////////////////////////\n");
	    // for(int i=0;i<wspace_size;i++)
	    // 	{
	    // 	  for(int j=0;j<=i;j++)
	    // 	    {
	    // 	      complex a={0,0};
	    // 	      for(int k=0;k<wspace_size;k++)
	    // 		for(int l=0;l<wspace_size;l++)
	    // 		  {
	    // 		    complex t;
	    // 		    unsafe_complex_conj1_prod(t,red_eig_vec[i+k*wspace_size],M[l+wspace_size*k]);
	    // 		    complex_summ_the_prod(a,t,red_eig_vec[j+l*wspace_size]);
	    // 		  }
	    // 	      master_printf("(%+4.4g,%+4.4g)\t",a[RE],a[IM]);
	    // 	    }
	    // 	  master_printf("\n");
	    // 	}
	    
	    //combine the basis vector to get the best eigenvectors approximations
	    wspace_size=wspace_min_size;
	    combine_basis_to_restart(wspace_min_size,wspace_max_size,red_eig_vec,V,mat_size);
	    
	    //set M to be diagonal
	    memset(M,0,sizeof(complex)*wspace_max_size*wspace_max_size);
	    for(int i=0;i<wspace_size;i++)
	      complex_put_to_real(M[i+i*wspace_max_size],red_eig_val[i]);
	    
	    // master_printf("//////////////////////////// matr //////////////////////////\n");
	    // for(int i=0;i<wspace_size;i++)
	    // 	{
	    // 	  for(int j=0;j<=i;j++)
	    // 	      master_printf("(%+4.4g,%+4.4g)\t",M[j+wspace_max_size*i][RE],M[j+wspace_max_size*i][IM]);
	    // 	  master_printf("\n");
	    // 	}
	    // for(int j=0;j<wspace_size;j++)
	    //   {
	    //     imp_mat(temp,V[j]);
	    //     for(int i=j;i<wspace_size;i++)
	    //       {
	    //         scalar_prod(M[j+wspace_max_size*i],V[i],temp,buffer,mat_size);
	    //         complex_conj(M[i+wspace_max_size*j],M[j+wspace_max_size*i]);
	    //       }
	    //   }
	    
	    // master_printf("//////////////////////////// true matr //////////////////////////\n");
	    // for(int i=0;i<wspace_size;i++)
	    // 	{
	    // 	  for(int j=0;j<=i;j++)
	    // 	      master_printf("(%+4.4g,%+4.4g)\t",M[j+wspace_max_size*i][RE],M[j+wspace_max_size*i][IM]);
	    // 	  master_printf("\n");
	    // 	}
	    // 	    crash("");
	  }
	
	//if not all eigenvectors have converged
	if(neig_conv!=neig)
	  {
	    //Solve the correction equation
	    vector_reset(V[wspace_size]);
	    
	    //find the target theta
	    double p_theta=
	      (residue_norm<eps_tr and residue_norm<red_eig_val[0] and residue_norm_old>residue_norm)?
	      red_eig_val[0]:
	      tau;
	    modified_GS(residue,eig_vec,buffer,neig_conv+1,mat_size);
	    
	    //projected matrix operator
	    auto proj_imp_mat=[p_theta,imp_mat,mat_size,buffer,eig_vec,neig_conv](complex *y,complex *x)
	      {
		//y=(A-theta)*x
		imp_mat(y,x);
		double_vector_summassign_double_vector_prod_double((double*)y,(double*)x,-p_theta,2*mat_size);
		
		//p=eig_vec^dagger*y
		for(int i=0;i<neig_conv+1;i++)
		  {
		    complex p;
		    scalar_prod(p,eig_vec[i],y,buffer,mat_size);
		    complex_vector_subtassign_complex_vector_prod_complex(y,eig_vec[i],p,mat_size);
		  }
	      };
	    
	    //solve the correction
	    complex *v=V[wspace_size];
	    double it_tol=pow(toldecay,-solvestep);
	    solvestep++;
	    cg_solve(v,proj_imp_mat,residue,mat_size,mat_size_to_allocate,it_tol*it_tol,linit_max,10000000);
	    
	    //orthonormalize v to eig_vec
	    modified_GS(v,eig_vec,buffer,neig_conv+1,mat_size);
	    const int max_cgs_it=5;
	    //orthogonalize v to V
	    double alpha=iterated_classical_GS(v,mat_size,wspace_size,V,buffer,max_cgs_it);
	    double_vector_prodassign_double((double*)v,1.0/alpha,2*mat_size);
	    
	    //update interaction matrix
	    imp_mat(temp,v);
	    for(int i=0;i<wspace_size+1;i++)
	      {
		scalar_prod(M[wspace_size+wspace_max_size*i],V[i],temp,buffer,mat_size);
		complex_conj(M[i+wspace_max_size*wspace_size],M[wspace_size+wspace_max_size*i]);
	      }
	    
	    wspace_size++;
	    
	    iter++;
	  }
      }
    while(iter<niter_max and neig_conv<neig);
    
    //free workspace
    for(int i=0;i<wspace_max_size;i++) nissa_free(V[i]);
    nissa_free(buffer);
    nissa_free(residue);
    nissa_free(temp);
  }
}

#endif
