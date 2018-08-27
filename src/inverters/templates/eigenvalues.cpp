#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "eigenvalues.hpp"

// #ifdef HAVE_EIGEN_DENSE
 #include "eigen3/Eigen/Dense"
// #endif

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
    void modified_GS(complex *v,complex **V,complex *buffer,int nvec,int vec_size)
    {
      for(int i=0;i<nvec;i++)
	{
	  complex s;
	  scalar_prod(s,V[i],v,buffer,vec_size);
	  complex_vector_subtassign_complex_vector_prod_complex(v,V[i],s,vec_size);
	}
    }
    
    //find eigenvalues of M and sort them according to |\lambda_i-\tau|
    void eigenvalues_of_hermatr_find_all_and_sort(double *lambda,complex *M,int size,double tau)
    {
      //structure to diagonalize
      using namespace Eigen;
      SelfAdjointEigenSolver<MatrixXcd> solver;
      
      //fill the matrix to be diagonalized
      MatrixXcd matr(size,size);
      for(int i=0;i<size;i++)
	for(int j=0;j<=i;j++)
	  matr(i,j)=std::complex<double>(M[j+size*i][RE],+M[j+size*i][IM]);
      
      //diagonalize
      solver.compute(matr);
      // std::cout<<"//////////////////////////// matr //////////////////////////"<<std::endl;
      // std::cout<<matr<<std::endl;
      std::cout<<"/////////////////////////// eigve //////////////////////////"<<std::endl;
      std::cout<<solver.eigenvectors()<<std::endl;
      std::cout<<"/////////////////////////// eigva //////////////////////////"<<std::endl;
      std::cout<<solver.eigenvalues()<<std::endl;
      std::cout<<"/////////////////////////////////////////////////////////////////"<<std::endl;
      
      //sort the eigenvalues and eigenvectors
      std::vector<std::tuple<double,double,int>> ei;
      for(int i=0;i<size;i++)
	{
	  double lambda=solver.eigenvalues()(i);
	  ei.push_back(std::make_tuple(fabs(lambda-tau),lambda,i));
	}
      std::sort(ei.begin(),ei.end());
      
      //fill output
      for(int ieig=0;ieig<size;ieig++)
	{
	  //fill eigvalue
	  using std::get;
	  lambda[ieig]=get<1>(ei[ieig]);
	  
	  //get index of what must be put in i
	  int ori=get<2>(ei[ieig]);
	  
	  //fill eigvec
	  for(int j=0;j<size;j++)
	    {
	      M[ieig+size*j][RE]=solver.eigenvectors()(j,ori).real();
	      M[ieig+size*j][IM]=solver.eigenvectors()(j,ori).imag();
	      
	      master_printf("%lg,%lg\t",M[ieig+size*j][RE],M[ieig+size*j][IM]);
	    }
	  master_printf("\n");
	  
	  //master_printf("%d   %lg %lg  %d\n",ieig,get<0>(ei[ieig]),get<1>(ei[ieig]),get<2>(ei[ieig]));
	}
    }
  }
}
