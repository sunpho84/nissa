#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "remez_algorithm.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../new_types/float256.h"
#include "../../routines/ios.h"
#include "../../routines/math.h"

//solve the linear system A*x=b
void linear_system_solve(float_256 *A,float_256 *x,float_256 *b,int n)
{
  int exch[n];
  
  //find the max of the row (Linf norm)
  for(int i=0;i<n;i++) 
    {
      exch[i]=i;
      float_256 row_norm;
      float_256_put_to_zero(row_norm);
      for(int j=0;j<n;j++) 
        {
          float_256 q;
          float_256_abs(q,A[i*n+j]);
          if(float_256_is_smaller(row_norm,q)) float_256_copy(row_norm,q);
        }
      if(row_norm==0) crash("num row norm");
      float_256_from_double(x[i],1);
      float_256_div(x[i],x[i],row_norm);
    }  
  
  //loop over the rows
  for(int k=0;k<n-1;k++) 
    {
      //find the pivot
      float_256 big;
      float_256_put_to_zero(big);
      int ipiv=0;
      for(int i=k;i<n;i++) 
        {
          int ip=exch[i];
          int ipk=n*ip+k;
          float_256 size;
          float_256_abs(size,A[ipk]);
          float_256_prodassign(size,x[ip]);
          
          if(float_256_is_greater(size,big)) 
            {
              float_256_copy(big,size);
              ipiv=i;
            }
        }
      if(big[0]==0) crash("null big");
        
      swap_ints(exch[ipiv],exch[k]);
      
      //pivotize
      float_256 pivot;
      float_256_copy(pivot,A[n*exch[k]+k]);
      for(int i=k+1;i<n;i++) 
        {
          float_256_div(A[n*exch[i]+k],A[n*exch[i]+k],pivot); //compute x to normalize other rows
          for(int j=k+1;j<n;j++) //subtract from the other rows
            float_256_subt_the_prod(A[n*exch[i]+j],A[n*exch[i]+k],A[n*exch[k]+j]);
        }
    }
  if(A[n*exch[n-1]+n-1][0]==0) crash("last element null");
  
  //build solution
  for(int i=0;i<n;i++)
    {
      float_256_copy(x[i],b[exch[i]]);
      for(int j=0;j<i;j++) float_256_subt_the_prod(x[i],A[n*exch[i]+j],x[j]);
    }
  
  //normalize with stored x
  for(int i=n-1;i>=0;i--) 
    {
      for(int j=i+1;j<n;j++) float_256_subt_the_prod(x[i],A[n*exch[i]+j],x[j]);
      float_256_div(x[i],x[i],A[n*exch[i]+i]);
    }
}

//find nmax_err_points+1 extrema of Chebyshev polynomial
void rat_approx_finder_t::find_cheb()
{
  double coeff=(maximum-minimum)/(exp(1)-1);
  for(int i=0;i<nmax_err_points;i++)
    {
      float_256_from_double(xmax[i],minimum+(exp(0.5*(1-cos((M_PI*i)/nzero_err_points)))-1)*coeff);
      float_256_from_double(zero[i],minimum+(exp(0.5*(1-cos(M_PI*(2*i+1)/(2*nzero_err_points))))-1)*coeff);
    }
}

//set the steps
void rat_approx_finder_t::set_step()
{
  float_256_subt(step[0],zero[0],minimum_256);
  for(int i=1;i<nzero_err_points;i++) float_256_subt(step[i],zero[i],zero[i-1]);
}

//perform an accomodation step
void rat_approx_finder_t::new_step()
{
  float_256 yy[nmax_err_points];

  float_256 eclose;
  float_256 farther;
  float_256_from_double(eclose,1.0e30);
  float_256_from_double(farther,0);
  
  //set left extrema
  float_256 zero0;
  float_256_from_double(zero0,minimum);

  for(int i=0;i<nmax_err_points;i++) 
    {
      //takes extrema
      float_256 zero1,xm;
      float_256_copy(zero1,zero[i]);
      if(i==nmax_err_points-1) float_256_from_double(zero1,maximum);
      float_256_copy(xm,xmax[i]);
      
      //check if we can move in one of the dirs
      float_256 xn,ym,yn,q;
      get_abs_err(ym,xm);
      float_256_copy(q,step[i]);
      float_256_summ(xn,xm,q);
      if(float_256_is_smaller(xn,zero0)||!float_256_is_smaller(xn,zero1)) 	// Cannot skip over adjacent boundaries
	{
	  float_256_uminus(q,q);
	  float_256_copy(xn,xm);
	  float_256_copy(yn,ym);
	}
      else 
	{
	  get_abs_err(yn,xn);
	  if(float_256_is_smaller(yn,ym)) 
	    {
	      float_256_uminus(q,q);
	      float_256_copy(xn,xm);
	      float_256_copy(yn,ym);
	    }
	}
      
      //move until reaching barrier or decreasing error
      int istep=0,quit=0;
      while(!quit&&!float_256_is_smaller(yn,ym))
       {
	 istep++;
	 if(istep>9) quit=1;
	 else
	   {
	     float_256 a;
	     float_256_copy(ym,yn);
	     float_256_copy(xm,xn);
	     float_256_summ(a,xm,q);
	     if(float_256_is_equal(a,xm)||!float_256_is_greater(a,zero0)||!float_256_is_smaller(a,zero1)) quit=1;
	     else
	       {
		 float_256_copy(xn,a);
		 get_abs_err(yn,xn);
	       }
	   }
       }
      
      //copy new position and max
      float_256_copy(xmax[i],xm);
      float_256_copy(yy[i],ym);

      //search extream of error
      if(float_256_is_greater(eclose,ym)) float_256_copy(eclose,ym);
      if(float_256_is_smaller(farther,ym)) float_256_copy(farther,ym);
      
      float_256_copy(zero0,zero1);
    }

  verbosity_lv3_master_printf("eclose: %16.16lg, farther: %16.16lg, spread: %16.16lg\n",eclose[0],farther[0],spread[0]);

  //constants
  float_256 half_256;float_256_from_double(half_256,0.5);
  float_256 quarter_256;float_256_from_double(quarter_256,0.25);
  float_256 mone;float_256_from_double(mone,-1);
  
  //decrease step size if error spread increased
  float_256 q,xm;
  float_256_subt(q,farther,eclose);
  //relative error spread
  if(eclose[0]!=0) float_256_div(q,q,eclose);
  //spread is increasing: decrease step size
  if(!float_256_is_smaller(q,spread)) float_256_prod(delta,delta,half_256);
  
  float_256_copy(spread,q);

  for(int i=0;i<nmax_err_points-1;i++) 
     {
       float_256_copy(q,yy[i+1]);
       if(q[0]!=0)
	 {
	   float_256_div(q,yy[i],q);
	   float_256_summassign(q,mone);
	 }
       else float_256_from_double(q,0.0625);
       if (float_256_is_greater(q,quarter_256)) float_256_from_double(q,0.25);
       
       float_256 d;
       float_256_subt(d,xmax[i+1],xmax[i]);
       float_256_prodassign(q,d);
       float_256_prod(step[i],q,delta);
     }
  float_256_copy(step[nmax_err_points-1],step[nmax_err_points-2]);
  
  //insert new locations for the zeros
  for(int i=0;i<nzero_err_points;i++)
    {	
      float_256_subt(xm,zero[i],step[i]);
      if(float_256_is_greater(xm,minimum_256))
	if(float_256_is_smaller(xm,maximum_256))
	  {
	    if(!float_256_is_greater(xm,xmax[i]))
	      {
		float_256_summ(xm,xmax[i],zero[i]);
		float_256_prodassign(xm,half_256);
	      }
	    if(!float_256_is_smaller(xm,xmax[i+1]))
	      {
		float_256_summ(xm,xmax[i+1],zero[i]);
		float_256_prodassign(xm,half_256);
	      }
	    float_256_copy(zero[i],xm);
	  }
    }
}

//set the linear system to be solved
void rat_approx_finder_t::set_linear_system(float_256 *matr,float_256 *vec) 
{
  for(int i=0;i<nzero_err_points;i++)
    {
      float_256 y;
      func_to_approx(y,zero[i]);
      
      float_256 z;
      float_256_from_double(z,1);
      int k=0;
      for(int j=0;j<=degree;j++)
        {
          float_256_copy(matr[i*nzero_err_points+k++],z); //x
          float_256_prod(z,z,zero[i]);
        }
      float_256_from_double(z,1);
      for(int j=0;j<degree;j++)
        {
	  float_256 yz; 
          float_256_prod(yz,y,z); //-xy
          float_256_uminus(matr[i*nzero_err_points+k++],yz);
          float_256_prod(z,z,zero[i]);
        }
      float_256_prod(vec[i],z,y); //x2y
    }
}

//compute separately the numerator and denominator of the approximation
void rat_approx_finder_t::compute_num_den_approx(float_256 yn,float_256 yd,float_256 x)
{
  float_256_copy(yn,coeff[degree]);
  for(int i=degree-1;i>=0;i--) 
    {
      float_256 c;
      float_256_prod(c,x,yn);
      float_256_summ(yn,c,coeff[i]);
    }
  float_256_summ(yd,x,coeff[2*degree]);
  for(int i=2*degree-1;i>degree;i--)
    {
      float_256 c;
      float_256_prod(c,x,yd);
      float_256_summ(yd,c,coeff[i]);
    }
}

//compute the approximation
void rat_approx_finder_t::compute_approx(float_256 y,float_256 x) 
{
  //numerator and denominator
  float_256 yn,yd;
  compute_num_den_approx(yn,yd,x);
  float_256_div(y,yn,yd);
}

//compute the error
void rat_approx_finder_t::get_err(float_256 err,float_256 x) 
{
  //compute exact fun and approx
  float_256 fun,app;
  func_to_approx(fun,x);
  compute_approx(app,x);
  
  //subtract
  float_256_subt(err,fun,app);

  //normalize
  if(fun[0]!=0) float_256_div(err,err,fun);
}

//compute absolute relative error
void rat_approx_finder_t::get_abs_err(float_256 err,float_256 x)
{
  //compute
  get_err(err,x);
  
  //takes absolute value
  float_256_abs(err,err);
}

//calculate the roots of the approximation
void rat_approx_finder_t::root_find(float_256 *roots,float_256 *poles,float_256 cons) 
{
  //define parameters
  const double upper=1,lower=-100000,tol=1.e-20;

  //find the numerator root
  float_256 poly[nmax_err_points];
  for(int i=0;i<=degree;i++) float_256_copy(poly[i],coeff[i]);
  
  for(int i=degree-1;i>=0;i--) 
     {
       root_find_Newton(roots[i],poly,i+1,lower,upper,tol);
       if(roots[i][0]==0) crash("Failure to converge on root %ld/%d",i+1,degree);
       float_256_div(poly[0],poly[0],roots[i]);
       float_256_uminus(poly[0],poly[0]);
       for(int j=1;j<=i;j++)
	 {
	   float_256_subt(poly[j],poly[j-1],poly[j]);
	   float_256_div(poly[j],poly[j],roots[i]);
	 }
     }
  
  //find the denominator roots
  float_256_from_double(poly[degree],1);
  for(int i=0;i<degree;i++) float_256_copy(poly[i],coeff[degree+1+i]);
  
  for(int i=degree-1;i>=0;i--) 
     {
       root_find_Newton(poles[i],poly,i+1,lower,upper,tol);
       if(poles[i][0]==0) crash("Failure to converge on pole %ld/%d",i+1,degree);
       float_256_div(poly[0],poly[0],poles[i]);
       float_256_uminus(poly[0],poly[0]);
       for(int j=1;j<=i;j++)
	 {
	   float_256_subt(poly[j],poly[j-1],poly[j]);
	   float_256_div(poly[j],poly[j],poles[i]);
	 }
     }
  
  float_256_copy(cons,coeff[degree]);
}

//evaluate the polynomial
void rat_approx_finder_t::poly_eval(float_256 f,float_256 x,float_256 *poly,int size) 
{
  float_256_copy(f,poly[size]);
  for(int i=size-1;i>=0;i--)
    {
      float_256_prodassign(f,x);
      float_256_summassign(f,poly[i]);
    }
}

//evaluate the derivative
void rat_approx_finder_t::poly_der(float_256 df,float_256 x,float_256 *poly,int size) 
{
  float_256 size_256;
  float_256_from_double(size_256,size);
  float_256_prod(df,size_256,poly[size]);
  for(int i=size-1;i>0;i--)
    {
      float_256 i_256;
      float_256_from_double(i_256,i);
      float_256_prodassign(df,x);
      float_256_summ_the_prod(df,i_256,poly[i]);
    }
}

//Newton's method to calculate roots
void rat_approx_finder_t::root_find_Newton(float_256 rtn,float_256 *poly,int size,double x1,double x2,double acc) 
{
  const int nmax_iter=10000;
  
  //start in the middle
  float_256_from_double(rtn,0.5*(x1+x2));
  
  //loop to find root
  int iter=0;
  float_256 dx;
  do
     {
       float_256 f,df;
       poly_eval(f,rtn,poly,size);
       poly_der(df,rtn,poly,size);
       float_256_div(dx,f,df);
       float_256_subtassign(rtn,dx);
       
       iter++;
     }
  while((iter<nmax_iter)&&(fabs(dx[0])>=acc));
  
  if(iter==nmax_iter) crash("Maximum number of iterations exceeded in Newton_root");
}

//decompose in a partial expansion
void get_partial_fraction_expansion(float_256 *res,float_256 *poles,float_256 *roots,float_256 cons,int n)
{
  float_256 numerator[n],denominator[n];
  for(int i=0;i<n;i++) float_256_copy(res[i],roots[i]);
  
  //construct the polynomials explicitly
  float_256_from_double(numerator[0],1);
  float_256_from_double(denominator[0],1);
  for(int i=1;i<n;i++)
    {
      float_256_put_to_zero(numerator[i]);
      float_256_put_to_zero(denominator[i]);
    }
  
  for(int j=0;j<n;j++)
    {
      for(int i=n-1;i>=0;i--)
	{
	  float_256_prodassign(numerator[i],res[j]);
	  float_256_uminus(numerator[i],numerator[i]);
	  float_256_prodassign(denominator[i],poles[j]);
	  float_256_uminus(denominator[i],denominator[i]);
	  
	  if(i>0)
	    {
	      float_256_summassign(numerator[i],numerator[i-1]);
	      float_256_summassign(denominator[i],denominator[i-1]);
	    }
	}
    }

  //convert to proper fraction form, because now is in the form 1+n/d
  for(int i=0;i<n;i++) float_256_subtassign(numerator[i],denominator[i]);

  //find the residues of the partial fraction expansion and absorb the coefficients
  for(int i=0;i<n;i++)
    {
      float_256_put_to_zero(res[i]);
      for(int j=n-1;j>=0;j--)
	{
	  float_256_prodassign(res[i],poles[i]);
	  float_256_summassign(res[i],numerator[j]);
	}
      for(int j=n-1;j>=0;j--)
	if(i!=j)
	  {
	    float_256 temp;
	    float_256_subt(temp,poles[i],poles[j]);
	    float_256_div(res[i],res[i],temp);
	  }
      float_256_prodassign(res[i],cons);
    }
  
  //res now holds the residues
  for(int i=0;i<n;i++) float_256_uminus(poles[i],poles[i]);

  //move the ordering of the poles from smallest to largest
  for(int j=0;j<n;j++)
    {
      int small=j;
      for(int i=j+1;i<n;i++)
	if(float_256_is_smaller(poles[i],poles[small])) small=i;
      
      if(small!=j)
	{
	  float_256_swap(poles[small],poles[j]);
	  float_256_swap(res[small],res[j]);
	}
    }
}

//generate the rational approximation
double rat_approx_finder_t::generate_approx(float_256 *weights,float_256 *poles,float_256 cons,double ext_minimum,double ext_maximum,
					    int ext_degree,int ext_num,int ext_den)
{
  //copy from out the degree and expo
  minimum=ext_minimum;
  maximum=ext_maximum;
  degree=ext_degree;
  num=ext_num;
  den=ext_den;

  //set degree depending coeffs
  nzero_err_points=2*degree+1;
  nmax_err_points=nzero_err_points+1;

  //alocate arrays
  step=nissa_malloc("step",nmax_err_points,float_256);
  coeff=nissa_malloc("coeff",nzero_err_points,float_256);
  zero=nissa_malloc("zero",nzero_err_points,float_256);
  xmax=nissa_malloc("xmax",nmax_err_points,float_256);

  //set the extrema of the interval
  float_256_from_double(minimum_256,minimum);
  float_256_from_double(maximum_256,maximum);
  
  //set delta, initial spread and tolerance
  float_256_from_double(spread,1.e37);
  float_256_from_double(delta,0.25);
  float_256_from_double(approx_tolerance,1e-15);
  
  //set the initial guess and set step
  find_cheb();
  set_step();

  //iterate up to convergence
  int iter=0;
  do
    {
      // 1) set up the system to be solved
      float_256 matr[nzero_err_points*nzero_err_points],vec[nzero_err_points];
      set_linear_system(matr,vec);
            
      // 2) solve the system
      linear_system_solve(matr,coeff,vec,nzero_err_points);
      
      // 3) find maxima and minima
      if(float_256_is_smaller(delta,approx_tolerance)) crash("delta too small, try increasing precision");
      if(float_256_is_greater(spread,approx_tolerance)) new_step();
      
      iter++;
    }
  while(float_256_is_greater(spread,approx_tolerance));

  verbosity_lv2_master_printf("Converged in %d iters\n",iter);
  
  //get err at max
  float_256 err;
  get_abs_err(err,xmax[0]);
  
  //compute the roots
  float_256 roots[degree];
  root_find(roots,poles,cons);
  
  //decompose
  get_partial_fraction_expansion(weights,poles,roots,cons,degree);
  
  for(int j=0;j<degree;j++)
    verbosity_lv2_master_printf("Residue = %lg, Pole = %lg\n",weights[j][0],poles[j][0]);
  verbosity_lv2_master_printf("Const: %lg\n",cons[0]);
  
  nissa_free(step);
  nissa_free(zero);
  nissa_free(coeff);
  nissa_free(xmax);
  
  //return the maximum error in the approximation
  return err[0];
}

