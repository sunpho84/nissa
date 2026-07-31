#ifndef _INCLUDE_REMEZ_ALGORITIHM_HPP
#define _INCLUDE_REMEZ_ALGORITIHM_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "new_types/high_prec.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void linear_system_solve(std::vector<float_high_prec_t>& A,
			   std::vector<float_high_prec_t>& x,
			   const std::vector<float_high_prec_t>& b,
			   const int& n);
  
  void get_partial_fraction_expansion(std::vector<float_high_prec_t>& res,
				      std::vector<float_high_prec_t>& poles,
				      const std::vector<float_high_prec_t>& roots,
				      const float_high_prec_t& cons,
				      const int& n);
  
  struct rat_approx_finder_t
  {
    //errors
    float_high_prec_t eclose;
    float_high_prec_t farther;
    
    //approximation parameters
    std::vector<float_high_prec_t> coeff;
    
    //parameters for solution of the system
    double minimum,maximum;
    int degree,nmax_err_points,nzero_err_points;
    
    //numerator and denominator of the power
    int num,den;
    
    //variables used to calculate the approximation
    std::vector<float_high_prec_t> zero;
    std::vector<float_high_prec_t> xmax;
    std::vector<float_high_prec_t> step;
    std::vector<float_high_prec_t> yy;
    std::vector<float_high_prec_t> poly;
    double delta,approx_tolerance;
    float_high_prec_t spread;
    
    //initial values of zero, maximal and steps
    void find_cheb();
    void set_step();
    
    //iter routine
    void set_linear_system(std::vector<float_high_prec_t>& matr,
			   std::vector<float_high_prec_t>& vec) const
    {
      for(int i=0;i<nzero_err_points;i++)
	{
	  float_high_prec_t y;
	  y=func_to_approx(zero[i]);
	  float_high_prec_t z=1;
	  int k=0;
	  for(int j=0;j<=degree;j++)
	    {
	      matr[i*nzero_err_points+(k++)]=z;
	      z*=zero[i];
	    }
	  z=1.0;
	  for(int j=0;j<degree;j++)
	    {
	      matr[i*nzero_err_points+(k++)]=-y*z;
	      z*=zero[i];
	    }
	  vec[i]=z*y;
	}
    }
    
    void new_step(int iter);
    
    /// Calculate the roots of the approximation
    void root_find(std::vector<float_high_prec_t>& roots,
		   std::vector<float_high_prec_t>& poles,
		   float_high_prec_t& cons) const
    {
      std::vector<float_high_prec_t> poly(nmax_err_points);
      
      //define parameters
      const double upper=1,lower=-100000,tol=1.e-20;
      
      //find the numerator root
      for(int i=0;i<=degree;i++) poly[i]=coeff[i];
      for(int i=degree-1;i>=0;i--)
	{
	roots[i]=root_find_Newton(poly,i+1,lower,upper,tol);
	if(roots[i]==0.0)
	  CRASH("Failure to converge on root %d/%d",i+1,degree);
	poly[0]/=-roots[i];
	for(int j=1;j<=i;j++) poly[j]=(poly[j-1]-poly[j])/roots[i];
	}
      
      //find the denominator roots
      poly[degree]=1;
      for(int i=0;i<degree;i++) poly[i]=coeff[degree+1+i];
      
      for(int i=degree-1;i>=0;i--)
      {
	poles[i]=root_find_Newton(poly,i+1,lower,upper,tol);
	if(poles[i]==0.0) CRASH("Failure to converge on pole %d/%d",i+1,degree);
	poly[0]/=-poles[i];
	for(int j=1;j<=i;j++) poly[j]=(poly[j-1]-poly[j])/poles[i];
      }
      
      cons=coeff[degree];
    }
    
    /// Evaluate a polynomial or its derivative
    float_high_prec_t poly_eval(const float_high_prec_t& x,
				const std::vector<float_high_prec_t>& poly,
				const int& size) const
    {
      float_high_prec_t f=poly[size];
      for(int i=size-1;i>=0;i--)
	f=f*x+poly[i];
      
      return f;
    }
    
    float_high_prec_t poly_eval(float_high_prec_t x,float_high_prec_t *poly,int size);
    
    float_high_prec_t poly_der(const float_high_prec_t& x,
			       const std::vector<float_high_prec_t>& poly,
			       const int& size) const
    {
      float_high_prec_t df=size*poly[size];
      for(int i=size-1;i>0;i--) df=df*x+i*poly[i];
      
      return df;
    }
    
    //Newton's method to calculate roots
    float_high_prec_t root_find_Newton(const std::vector<float_high_prec_t>& poly,
				       const int& size,
				       const double& x1,
				       const double& x2,
				       const double& acc) const
    {
      const int nmax_iter=10000;
      
      //start in the middle
      float_high_prec_t rtn=0.5*(x1+x2);
      
      //loop to find root
      int iter=0;
      float_high_prec_t dx;
      do
	{
	  const float_high_prec_t f=poly_eval(rtn,poly,size);
	  const float_high_prec_t df=poly_der(rtn,poly,size);
	  dx=f/df;
	  
	  rtn-=dx;
	  
	  iter++;
	}
      while((iter<nmax_iter)&&(abs(dx)>=acc));
      
      if(iter==nmax_iter) CRASH("Maximum number of iterations exceeded in Newton_root");
      
      return rtn;
    }
    
    //calculate function required for the approximation, or its error
    float_high_prec_t func_to_approx(const float_high_prec_t& x) const
    {
      return float_high_prec_t_pow_int_frac(x,num,den);
    }
    
    float_high_prec_t get_abs_err(const float_high_prec_t& x) const
    {
      return abs(get_err(x));
    }
    
    float_high_prec_t get_err(float_high_prec_t x) const
    {
      //compute exact fun and approx
      const float_high_prec_t fun=func_to_approx(x);
      const float_high_prec_t app=compute_approx(x);
      
      //subtract
      float_high_prec_t err=
	fun-app;
      
      //normalize
      if(fun.get_d()!=0)
	err/=fun;
      
      return err;
    }
    
    /// Compute separately the numerator and denominator of the approximation
    void compute_num_den_approx(float_high_prec_t &yn,
				float_high_prec_t &yd,
				const float_high_prec_t& x) const
    {
      yn=coeff[degree];
      for(int i=degree-1;i>=0;i--)
	yn=x*yn+coeff[i];
      yd=x+coeff[2*degree];
      for(int i=2*degree-1;i>degree;i--)
	yd=x*yd+coeff[i];
    }
    
    float_high_prec_t compute_approx(const float_high_prec_t& x) const
    {
      //numerator and denominator
      float_high_prec_t yn,yd;
      compute_num_den_approx(yn,yd,x);
      return yn/yd;
    }
    
    //generate the rational approximation with a given number of poles,
    //checking that min and max err are within a factor of 1+toll, and
    //giving up if min comes out to be larger than target_err
    double generate_approx(std::vector<float_high_prec_t>& weights,
			   std::vector<float_high_prec_t>& poles,
			   float_high_prec_t &cons,
			   const double& ext_minimum,
			   const double& ext_maximum,
			   const int& ext_degree,
			   const int& ext_num,
			   const int& ext_den,
			   const double& target_err,
			   const double& toll)
    {
      
      //if target_err is not positive, it is ignored
      const bool consider_err=(target_err>0);
      
      //copy from out the degree and expo
      minimum=ext_minimum;
      maximum=ext_maximum;
      degree=ext_degree;
      num=ext_num;
      den=ext_den;
      
      //set degree depending coeffs
      nzero_err_points=2*degree+1;
      nmax_err_points=nzero_err_points+1;
      
      //set delta, initial spread and tolerance
      spread=1e37;
      delta=0.25;
      approx_tolerance=toll;
      
      //alocate arrays
      std::vector<float_high_prec_t> matr(nzero_err_points*nzero_err_points);
      std::vector<float_high_prec_t> vec(nzero_err_points);
      step.resize(nmax_err_points);
      coeff.resize(nmax_err_points);
      zero.resize(nmax_err_points);
      xmax.resize(nmax_err_points);
      
      //set the initial guess and set step
      find_cheb();
      set_step();
      
      //iterate up to convergence
      int iter=0;
      do
	{
	  // 1) set up the system to be solved
	  set_linear_system(matr,vec);
	  
	  // 2) solve the system
	  linear_system_solve(matr,coeff,vec,nzero_err_points);
	  
	  // 3) find maxima and minima
	  if(iter==0 || (spread>approx_tolerance && (farther>target_err || not consider_err))) new_step(iter);
	  
	  if(delta<approx_tolerance)
	    {
	      WARNING("reached precision %lg while computing %d terms approximation of x^(%d/%d) with tolerance %lg\n",
		      spread.get_d(),degree,num,den,approx_tolerance);
	      MASTER_PRINTF("precision not enough to reach %lg precision requested!!!\n",approx_tolerance);
#if HIGH_PREC_TYPE==NATIVE_HIGH_PREC
	      MASTER_PRINTF("use GMP if possible!\n");
#else
	      MASTER_PRINTF("compile with higher precision!\n");
#endif
	    }
	  iter++;
	}
      //while(float_high_prec_t_is_greater(spread,approx_tolerance));
      while(spread>approx_tolerance && delta>=approx_tolerance && ((not consider_err) || (eclose.get_d()<=target_err && farther.get_d()>target_err)));
      
      //write some info
      if(spread<=approx_tolerance) VERBOSITY_LV3_MASTER_PRINTF("Spread %lg reduced below %lg\n",spread.get_d(),approx_tolerance);
      if(consider_err && eclose>target_err)  VERBOSITY_LV3_MASTER_PRINTF("Accuracy cannot be better than %lg when %lg asked\n",eclose.get_d(),target_err);
      
      //get err at max and check
      if((consider_err and farther.get_d()<=target_err) or ((not consider_err) and (spread<=approx_tolerance)))
	{
	  VERBOSITY_LV2_MASTER_PRINTF("Converged with %d zeroes in %d iters, maxerr %lg when asked %lg\n",nzero_err_points,iter,farther.get_d(),target_err);
	  
	  //compute the roots
	  std::vector<float_high_prec_t> roots(degree);
	  root_find(roots,poles,cons);
	  
	  //decompose
	  get_partial_fraction_expansion(weights,poles,roots,cons,degree);
	  
	  for(int j=0;j<degree;j++)
	    VERBOSITY_LV2_MASTER_PRINTF("Weight = %16.16lg, Pole = %16.16lg\n",weights[j].get_d(),poles[j].get_d());
	  VERBOSITY_LV2_MASTER_PRINTF("Const: %16.16lg\n",cons.get_d());
	}
      else VERBOSITY_LV2_MASTER_PRINTF("Not converged to %lg prec with %d poles in %d iters (reached: %lg)\n",target_err,degree,iter,farther.get_d());
      
      //return the maximum error in the approximation
      const double ret=
	farther.get_d();
      
      return ret;
    }
  };
  
  double generate_approx(rat_approx_t& appr,
			 const double& minimum,
			 const double& maximum,
			 const int& num,
			 const int& den,
			 const double& minerr,
			 const double& tollerance);
  
  void generate_approx_of_maxerr(rat_approx_t &appr,
				 const double& minimum,
				 const double& maximum,
				 const double& maxerr,
				 const int& num,
				 const int& den,
				 const char *name="");
}

#endif
