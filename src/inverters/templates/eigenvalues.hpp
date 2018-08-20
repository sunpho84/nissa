#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

namespace nissa
{
  //class to get eigenvalues from a matrix op
  template <class T,class F1,class F2>
  class eigenvalues_herm_finder_t
  {
    //matrix size
    int mat_size;
    
    //stopping criterion
    double tol;
    
    //minimum or maximum
    bool min_max;
    
    //function to fill the random source
    const F1 &fill;
    
    //operator implementing the hermitian matrix to be solved
    const F2 &op;
    
    //max number ofiterations
    int nmax_iter{100000};
    
  public:
    
    //constructor
    eigenvalues_herm_finder_t<T,F1,F2>(int mat_size,double tol,bool min_max,const F1 &fill,const F2 &op) : mat_size(mat_size),tol(tol),min_max(min_max),fill(fill),op(op)
    {
      
    }
    
    //destructor
    ~eigenvalues_herm_finder_t()
    {
    }
    
    //set the maximal number of iterations
    void set_nmax_iter(int n)
    {
      nmax_iter=n;
    }
    
    //compute eigenvectors
    void get(T **eig_vec,double *eig_val,int neig)
    {
    }
  };
  
  //helper function to avoid explicitating the types
  template <class T,class F1,class F2>
  auto get_eigenvalues_herm_finder(int mat_size,double tol,bool min_max,const F1 &fill,const F2 &op)
  {
    return eigenvalues_herm_finder_t<T,decltype(fill),decltype(op)>(mat_size,tol,min_max,fill,op);
  }
}

#endif
