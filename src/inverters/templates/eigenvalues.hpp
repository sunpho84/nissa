#ifndef _EIGENVALUES_HPP
#define _EIGENVALUES_HPP

namespace nissa
{
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
    
  public:
    //constructor
    eigenvalues_herm_finder_t<T,F1,F2>(int mat_size,double tol,bool min_max,const F1 &fill,const F2 &op) : mat_size(mat_size),tol(tol),min_max(min_max),fill(fill),op(op)
    {
      
    }
    
    //compute eigenvectors
    void get(T **eig_vec,double *eig_val,int neig)
    {
    }
  };
}

#endif
