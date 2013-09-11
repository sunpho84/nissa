#ifndef _INCLUDE_REMEZ_ALGORITIHM_H
#define _INCLUDE_REMEZ_ALGORITIHM_H

#include "new_types/new_types_definitions.h"
#include "new_types/float256.h"
#include "base/vectors.h"

class rat_approx_finder_t
{
 private:
  //approximation parameters
  float_256 *coeff;

  //parameters for solution of the system
  double minimum,maximum;
  float_256 minimum_256,maximum_256;
  int degree,nmax_err_points,nzero_err_points;
  
  //numerator and denominator of the power
  int num,den;

  //variables used to calculate the approximation
  float_256 *zero,*xmax,*step;
  float_256 delta,spread,approx_tolerance;

  //initial values of zero, maximal and steps
  void find_cheb();
  void set_step();
  
  //iter routine
  void set_linear_system(float_256 *matr,float_256 *vec);
  void new_step(); 
  
  //calculate the roots of the approximation
  void root_find(float_256 *roots,float_256 *poles,float_256 cons);

  //evaluate a polynomial or its derivative
  void poly_eval(float_256 out,float_256 x,float_256 *poly,int size);
  void poly_der(float_256 out,float_256 x,float_256 *poly,int size);

  //Newton's method to calculate roots
  void root_find_Newton(float_256 out,float_256 *poly,int i,double x1,double x2,double acc);

  //calculate function required for the approximation, or its error
  void func_to_approx(float_256 y,float_256 x){float_256_pow_int_frac(y,x,num,den);}
  void get_abs_err(float_256 err,float_256 x);
  void get_err(float_256 err,float_256 x);

  //compute the approximation
  void compute_num_den_approx(float_256 yn,float_256 yd,float_256 x);
  void compute_approx(float_256 out,float_256 x);

 public:  
  //generate the rational approximation
  double generate_approx(float_256 *weights,float_256 *pole,float_256 cons,double ext_minimum,double ext_maximum,int ext_degree,int num,int den);
};


#endif


