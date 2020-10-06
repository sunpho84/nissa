#ifndef _INCLUDE_REMEZ_ALGORITIHM_HPP
#define _INCLUDE_REMEZ_ALGORITIHM_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/high_prec.hpp"

#include "base/vectors.hpp"

namespace nissa
{
  float_high_prec_t float_high_prec_t_pow_int_frac(float_high_prec_t ext_in,int n,int d);
  
  struct rat_approx_finder_t
  {
    //errors
    float_high_prec_t eclose;
    float_high_prec_t farther;
    
    //approximation parameters
    float_high_prec_t *coeff;
    
    //parameters for solution of the system
    double minimum,maximum;
    int degree,nmax_err_points,nzero_err_points;
    
    //numerator and denominator of the power
    int num,den;
    
    //variables used to calculate the approximation
    float_high_prec_t *zero,*xmax,*step,*yy,*poly;
    double delta,approx_tolerance;
    float_high_prec_t spread;
    
    //initial values of zero, maximal and steps
    void find_cheb();
    void set_step();
    
    //iter routine
    void set_linear_system(float_high_prec_t *matr,float_high_prec_t *vec);
    void new_step(int iter);
    
    //calculate the roots of the approximation
    void root_find(float_high_prec_t *roots,float_high_prec_t *poles,float_high_prec_t &cons);
    
    //evaluate a polynomial or its derivative
    float_high_prec_t poly_eval(float_high_prec_t x,float_high_prec_t *poly,int size);
    float_high_prec_t poly_der(float_high_prec_t x,float_high_prec_t *poly,int size);
    
    //Newton's method to calculate roots
    float_high_prec_t root_find_Newton(float_high_prec_t *poly,int i,double x1,double x2,double acc);
    
    //calculate function required for the approximation, or its error
    float_high_prec_t func_to_approx(float_high_prec_t x)
    {return float_high_prec_t_pow_int_frac(x,num,den);}
    float_high_prec_t get_abs_err(float_high_prec_t x);
    float_high_prec_t get_err(float_high_prec_t x);
    
    //compute the approximation
    void compute_num_den_approx(float_high_prec_t &yn,float_high_prec_t &yd,float_high_prec_t x);
    float_high_prec_t compute_approx(float_high_prec_t x);
    
    //generate the rational approximation
    double generate_approx(float_high_prec_t *weights,float_high_prec_t *pole,float_high_prec_t &cons,double ext_minimum,double ext_maximum,int ext_degree,int num,int den,double minerr,double tollerance);
  };
  double generate_approx(rat_approx_t &appr,double minimum,double maximum,int num,int den,double minerr,double tollerance);
  void generate_approx_of_maxerr(rat_approx_t &appr,double minimum,double maximum,double maxerr,int num,int den,const char *name="");
}

#endif
