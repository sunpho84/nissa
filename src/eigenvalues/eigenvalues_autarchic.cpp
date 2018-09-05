#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "eigenvalues_autarchic.hpp"

#ifdef USE_EIGEN
 #include <eigen3/Eigen/Dense>
#endif

//part of this should be moved to linalgs?

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
    
    //Ai=Ai+Bi*c
    void complex_vector_summassign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n)
    {
      GET_THREAD_ID();
      
      NISSA_PARALLEL_LOOP(i,0,n)
	complex_summ_the_prod(a[i],b[i],c);
      set_borders_invalid(a);
    }
    
    //orthogonalize v with respect to the nvec of V
    void modified_GS(complex *v,complex **V,complex *buffer,int nvec,int vec_length)
    {
      for(int i=0;i<nvec;i++)
	{
	  complex s;
	  scalar_prod(s,V[i],v,buffer,vec_length);
	  complex_vector_subtassign_complex_vector_prod_complex(v,V[i],s,vec_length);
	}
    }
    
    //find eigenvalues of M and sort them according to |\lambda_i-\tau|
    //NB: M is larger than neig
    void eigenvalues_of_hermatr_find_all_and_sort(complex *eig_vec,double *lambda,const complex *M,const int M_size,const int neig,const double tau)
    {
#if !USE_EIGEN
      crash("need Eigen");
#else
      //structure to diagonalize
      using namespace Eigen;
      SelfAdjointEigenSolver<MatrixXcd> solver;
      
      //fill the matrix to be diagonalized
      // master_printf("//////////////////////////// matr //////////////////////////\n");
      MatrixXcd matr(neig,neig);
      for(int i=0;i<neig;i++)
	for(int j=0;j<=i;j++)
	  matr(i,j)=std::complex<double>(M[j+M_size*i][RE],M[j+M_size*i][IM]);
      
      //diagonalize
      solver.compute(matr);
      
      //sort the eigenvalues and eigenvectors
      std::vector<std::tuple<double,double,int>> ei;
      for(int i=0;i<neig;i++)
	{
	  double lambda=solver.eigenvalues()(i);
	  ei.push_back(std::make_tuple(fabs(lambda-tau),lambda,i));
	}
      std::sort(ei.begin(),ei.end());
      
      //fill output
      for(int ieig=0;ieig<neig;ieig++)
	{
	  //fill eigvalue
	  using std::get;
	  lambda[ieig]=get<1>(ei[ieig]);
	  
	  //get index of what must be put in i
	  int ori=get<2>(ei[ieig]);
	  
	  //fill eigvec
	  for(int j=0;j<neig;j++)
	    {
	      eig_vec[ieig+neig*j][RE]=solver.eigenvectors()(j,ori).real();
	      eig_vec[ieig+neig*j][IM]=solver.eigenvectors()(j,ori).imag();
	    }
	}
#endif
    }
    
    double iterated_classical_GS(complex *v,int vec_length,int nvec,complex **A,complex *buffer,const int max_cgs_it)
    {
      const double alpha=0.5;
      bool isorth=0;
      double v_norm=0.0;
      double v_norm_old=sqrt(double_vector_glb_norm2(v,vec_length));
      
      for(int i=0;i<max_cgs_it and not isorth;i++)
	{
	  for(int j=0;j<nvec;j++)
	    {
	      complex p;
	      scalar_prod(p,A[j],v,buffer,vec_length);
	      complex_vector_subtassign_complex_vector_prod_complex(v,A[j],p,vec_length);
	    }
	  v_norm=sqrt(double_vector_glb_norm2(v,vec_length));
	  
	  isorth=(v_norm>alpha*v_norm_old);
	  v_norm_old=v_norm;
	}
      
      return v_norm;
    }
    
    //form v[0:nout]=v[0:nin]*coeffs[0:nin,0:nout], using v itself
    void combine_basis_to_restart(int nout,int nin,complex *coeffs,complex **vect,int vec_length)
    {
      GET_THREAD_ID();
      
      NISSA_PARALLEL_LOOP(iel,0,vec_length)
	{
	  //store the linear combo at fixed i
	  complex tmp[nout];
	  
	  //combine
	  for(int j=0;j<nout;j++)
	    {
	      complex_put_to_zero(tmp[j]);
	      for(int k=0;k<nin;k++)
		complex_summ_the_prod(tmp[j],vect[k][iel],coeffs[j+nin*k]);
	    }
	  
	  //overwrite
	  for(int j=0;j<nout;j++)
	    complex_copy(vect[j][iel],tmp[j]);
	}
      
      //barrier and mark as modified
      for(int j=0;j<nout;j++)
	set_borders_invalid(vect[j]);
    }
  }
}
