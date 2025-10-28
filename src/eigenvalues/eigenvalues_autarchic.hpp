#ifndef _EIGENVALUES_AUTARCHIC_HPP
#define _EIGENVALUES_AUTARCHIC_HPP

#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "inverters/templates/modern_cg.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  namespace internal_eigenvalues
  {
    void modified_GS(complex *v,complex **V,int nvec,int vec_size);
    double iterated_classical_GS(complex *v,int vec_size,int nvec,complex **A,const int max_cgs_it);
    void eigenvalues_find_all_and_sort(complex *eig_vec,int eig_vec_row_size,double *lambda,const complex *M,const int M_size,const int neig,const double tau);
    void combine_basis_to_restart(int nout,int nin,complex *coeffs,int coeffs_row_length,complex **vect,int vec_length);
  }
  
  //reimplementation of the adaptation made by Carsten Urbach of Jacobi-Davidson algorithm by R. Geus and O. Chinellato
  
  //find the neig eigenvalues closest to the target
  template <class Fmat,class Filler>
  void eigenvalues_find_autarchic(complex **eig_vec,complex *eig_val,int neig,bool min_max,
					     const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
					     const double target_precision,const int niter_max,
					     const Filler &filler)
  {
    CRASH("reimplement");

    // using namespace internal_eigenvalues;
    
    // if(neig>mat_size) CRASH("Asking to find %d eigenvectors for a matrix of rank %d",neig,mat_size);
    
    // //set pars
    // const double tau_list[2]={0.0,50.0};
    // double tau=tau_list[min_max];
    // const int linit_max=200;
    // const double toldecay_list[2]={1.7,1.5};
    // const double toldecay=toldecay_list[min_max];
    // const double eps_tr_list[2]={1e-3,5e-2};
    // const double eps_tr=eps_tr_list[min_max];
    // const int wspace_min_size=std::max(8,neig);
    // const int wspace_max_size=2*wspace_min_size;
    
    // //allocate workspace
    // complex *V[wspace_max_size];
    // for(int i=0;i<wspace_max_size;i++) V[i]=nissa_malloc("Vi",mat_size_to_allocate,complex);
    // complex *residue=nissa_malloc("residue",mat_size_to_allocate,complex);
    // complex *temp=nissa_malloc("temp",mat_size,complex);
    // double *red_eig_val=new double[wspace_max_size];
    // complex *red_eig_vec=new complex[wspace_max_size*wspace_max_size];
    
    // //fill V and orthonormalize
    // for(int i=0;i<wspace_max_size;i++)
    //   {
    // 	//fill
    // 	filler(V[i]);
    // 	//orthogonalize
    // 	modified_GS(V[i],V,i,mat_size);
    // 	//normalize
    // 	double useless_rat;
    // 	double_vector_normalize(&useless_rat,(double*)(V[i]),(double*)(V[i]),1.0,2*mat_size);
    //   }
    
    // //generate interaction matrix
    // int wspace_size=wspace_min_size;
    // complex *M=new complex[wspace_max_size*wspace_max_size];
    // for(int j=0;j<wspace_size;j++)
    //   {
    // 	imp_mat(temp,V[j]);
    // 	for(int i=j;i<wspace_size;i++)
    // 	  {
    // 	    complex_vector_glb_scalar_prod(M[j+wspace_max_size*i],V[i],temp,mat_size);
    // 	    complex_conj(M[i+wspace_max_size*j],M[j+wspace_max_size*i]);
    // 	    VERBOSITY_LV3_MASTER_PRINTF("M_(%d,%d)=%lg,%lg\n",j,i,M[i+wspace_max_size*j][RE],M[i+wspace_max_size*j][IM]);
    // 	  }
    //   }
    
    // //main loop
    // int solvestep=1;
    // int iter=0;
    // int neig_conv=0;
    // do
    //   {
    // 	MASTER_PRINTF("Iteration %d, wspace size: %d [%d:%d]\n",iter,wspace_size,wspace_min_size,wspace_max_size);
	
    // 	//reset residue norm
    // 	double residue_norm=0.0;
	
    // 	//find all eigenvalues of the reduced problem, sort them by distance with tau
    // 	eigenvalues_find_all_and_sort(red_eig_vec,wspace_max_size,red_eig_val,M,wspace_max_size,wspace_size,tau);
	
    // 	//combine the vectors
    // 	complex *e=eig_vec[neig_conv];
    // 	vector_reset(e);
    // 	for(int j=0;j<wspace_size;j++)
    // 	  complex_vector_summassign_complex_vector_prod_complex(e,V[j],red_eig_vec[wspace_max_size*j],mat_size);
	
    // 	//compute the residue
    // 	double residue_norm_old=residue_norm;
    // 	imp_mat(residue,e);
    // 	double_vector_summassign_double_vector_prod_double((double*)residue,(double*)e,-red_eig_val[0],mat_size*2);
    // 	residue_norm=sqrt(double_vector_glb_norm2(residue,mat_size));
    // 	MASTER_PRINTF("eig: %.16lg, res: %.16lg\n",red_eig_val[0],residue_norm);
	
    // 	//if converged
    // 	bool this_converged=(residue_norm<target_precision);
	
    // 	if(this_converged)
    // 	  {
    // 	    //store eigenvalue
    // 	    complex_put_to_real(eig_val[neig_conv],red_eig_val[0]);
    // 	    MASTER_PRINTF("Eigenvalue %d/%d, %lg converged!\n",neig_conv,neig,eig_val[neig_conv][RE]);
	    
    // 	    //shift the others back
    // 	    for(int i=0;i<wspace_size-1;i++)
    // 	      red_eig_val[i]=red_eig_val[i+1];
	    
    // 	    //restart using the remaining eigenvectors as basis
    // 	    complex *coeffs=red_eig_vec+1; //Set a shifted view on red_eig_vec
    // 	    combine_basis_to_restart(wspace_size-1,wspace_size,coeffs,wspace_max_size,V,mat_size);
	    
    // 	    //update workspace size
    // 	    wspace_size--;
	    
    // 	    //make M diagonal
    // 	    for(int i=0;i<wspace_size;i++)
    // 	      for(int j=0;j<wspace_size;j++)
    // 		{
    // 		  const double val=(i==j)?red_eig_val[i]:0.0;
    // 		  complex_put_to_real(M[j+i*wspace_max_size],val);
    // 		}
	    
    // 	    //set tau closer to the minimal eig_val
    // 	    if(min_max==0)
    // 	      tau=std::max(eig_val[neig_conv][RE],tau);
    // 	    else
    // 	      tau=std::min(eig_val[neig_conv][RE],tau);
	    
    // 	    //reset the stopping criterion
    // 	    solvestep=1;
	    
    // 	    //number converged incremented at the end of the loop
    // 	  }
	
    // 	//reset if exceeded the workspace size
    // 	if(wspace_size==wspace_max_size)
    // 	  {
    // 	    MASTER_PRINTF("Resetting\n");
	    
    // 	    //combine the basis vector to get the best eigenvectors approximations
    // 	    wspace_size=wspace_min_size;
    // 	    combine_basis_to_restart(wspace_min_size,wspace_max_size,red_eig_vec,wspace_max_size,V,mat_size);
	    
    // 	    //set M to be diagonal
    // 	    for(int i=0;i<wspace_size;i++)
    // 	      for(int j=0;j<wspace_size;j++)
    // 		{
    // 		  const double val=(i==j)?red_eig_val[i]:0.0;
    // 		  complex_put_to_real(M[j+i*wspace_max_size],val);
    // 		}
    // 	  }
	
    // 	//if not all eigenvectors have converged
    // 	if(neig_conv!=neig)
    // 	  {
    // 	    //Solve the correction equation
    // 	    vector_reset(V[wspace_size]);
	    
    // 	    //find the target theta
    // 	    double p_theta=
    // 	      (residue_norm<eps_tr and residue_norm<red_eig_val[0] and residue_norm_old>residue_norm)?
    // 	      red_eig_val[0]:
    // 	      tau;
    // 	    modified_GS(residue,eig_vec,neig_conv+1,mat_size);
	    
    // 	    //projected matrix operator
    // 	    auto proj_imp_mat=[p_theta,imp_mat,mat_size,eig_vec,neig_conv](complex *y,complex *x)
    // 	      {
    // 		//y=(A-theta)*x
    // 		imp_mat(y,x);
    // 		double_vector_summassign_double_vector_prod_double((double*)y,(double*)x,-p_theta,2*mat_size);
		
    // 		//p=eig_vec^dagger*y
    // 		for(int i=0;i<neig_conv+1;i++)
    // 		  {
    // 		    complex p;
    // 		    complex_vector_glb_scalar_prod(p,eig_vec[i],y,mat_size);
    // 		    complex_vector_subtassign_complex_vector_prod_complex(y,eig_vec[i],p,mat_size);
    // 		  }
    // 	      };
	    
    // 	    //solve the correction
    // 	    complex *v=V[wspace_size];
    // 	    double it_tol=pow(toldecay,-solvestep);
    // 	    solvestep++;
    // 	    cg_solve(v,proj_imp_mat,residue,mat_size,mat_size_to_allocate,it_tol*it_tol,linit_max,10000000);
	    
    // 	    //orthonormalize v to eig_vec
    // 	    modified_GS(v,eig_vec,neig_conv+1,mat_size);
    // 	    const int max_cgs_it=5;
    // 	    //orthogonalize v to V
    // 	    double alpha=iterated_classical_GS(v,mat_size,wspace_size,V,max_cgs_it);
    // 	    double_vector_prodassign_double((double*)v,1.0/alpha,2*mat_size);
    // 	    MASTER_PRINTF("alpha: %.16lg\n",alpha);
	    
    // 	    //update interaction matrix
    // 	    imp_mat(temp,v);
    // 	    for(int i=0;i<wspace_size+1;i++)
    // 	      {
    // 		complex_vector_glb_scalar_prod(M[wspace_size+wspace_max_size*i],V[i],temp,mat_size);
    // 		complex_conj(M[i+wspace_max_size*wspace_size],M[wspace_size+wspace_max_size*i]);
    // 	      }
	    
    // 	    wspace_size++;
	    
    // 	    iter++;
	    
    // 	    //increment only here not to scramble the counting
    // 	    if(this_converged)
    // 	      neig_conv++;
    // 	  }
    //   }
    // while(iter<niter_max and neig_conv<neig);
    
    // //free workspace
    // delete [] red_eig_val;
    // delete [] red_eig_vec;
    // delete [] M;
    // for(int i=0;i<wspace_max_size;i++) nissa_free(V[i]);
    // nissa_free(residue);
    // nissa_free(temp);
  }
}

#endif
