#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "eigenvalues_autarchic.hpp"

#ifdef USE_EIGEN
 #include <Eigen/Dense>
#endif

namespace nissa
{
  namespace internal_eigenvalues
  {
    //orthogonalize v with respect to the nvec of V
    void modified_GS(complex *v,complex **V,int nvec,int vec_length)
    {
      for(int i=0;i<nvec;i++)
	{
	  complex s;
	  complex_vector_glb_scalar_prod(s,V[i],v,vec_length);
	  complex_vector_subtassign_complex_vector_prod_complex(v,V[i],s,vec_length);
	}
    }
    
    //find eigenvalues of M and sort them according to |\lambda_i-\tau|
    //NB: M is larger than neig
    void eigenvalues_find_all_and_sort(complex *eig_vec,int eig_vec_row_size,double *lambda,const complex *M,const int M_size,const int neig,const double tau)
    {
#if !USE_EIGEN
      CRASH("need Eigen");
#else
      //structure to diagonalize
      using namespace Eigen;
      SelfAdjointEigenSolver<MatrixXcd> solver;
      
      //fill the matrix to be diagonalized
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
	      eig_vec[ieig+eig_vec_row_size*j][RE]=solver.eigenvectors()(j,ori).real();
	      eig_vec[ieig+eig_vec_row_size*j][IM]=solver.eigenvectors()(j,ori).imag();
	    }
	}
#endif
    }
    
    double iterated_classical_GS(complex *v,int vec_length,int nvec,complex **A,const int max_cgs_it)
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
	      complex_vector_glb_scalar_prod(p,A[j],v,vec_length);
	      complex_vector_subtassign_complex_vector_prod_complex(v,A[j],p,vec_length);
	    }
	  v_norm=sqrt(double_vector_glb_norm2(v,vec_length));
	  
	  isorth=(v_norm>alpha*v_norm_old);
	  v_norm_old=v_norm;
	}
      
      return v_norm;
    }
    
    //form v[0:nout]=v[0:nin]*coeffs[0:nin,0:nout], using v itself
    void combine_basis_to_restart(int nout,int nin,complex *coeffs,int coeffs_row_length,complex **vect,int vec_length)
    {
      
      NISSA_PARALLEL_LOOP(iel,0,vec_length)
	{
	  //store the linear combo at fixed i
	  complex *tmp=new complex[nout];
	  
	  //combine
	  for(int j=0;j<nout;j++)
	    {
	      complex_put_to_zero(tmp[j]);
	      for(int k=0;k<nin;k++)
		complex_summ_the_prod(tmp[j],vect[k][iel],coeffs[j+coeffs_row_length*k]);
	    }
	  
	  //overwrite
	  for(int j=0;j<nout;j++)
	    complex_copy(vect[j][iel],tmp[j]);
	  
	  delete[] tmp;
	}
      NISSA_PARALLEL_LOOP_END;
      
      //barrier and mark as modified
      for(int j=0;j<nout;j++)
	set_borders_invalid(vect[j]);
    }
  }
}
