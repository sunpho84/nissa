#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "new_types/high_prec.hpp"
#include "new_types/rat_approx.hpp"
#include "remez_algorithm.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //solve the linear system A*x=b
  void linear_system_solve(float_high_prec_t *A,float_high_prec_t *x,float_high_prec_t *b,int n)
  {
    int exch[n];
    
    //find the max of the row (Linf norm)
    for(int i=0;i<n;i++) 
      {
	exch[i]=i;
	float_high_prec_t row_norm=0.0;
	for(int j=0;j<n;j++) 
	  {
	    float_high_prec_t q;
	    q=abs(A[i*n+j]);
	    if(row_norm<q) row_norm=q;
	  }
	if(row_norm==0.0) crash("num row norm");
	x[i]=1/row_norm;
      }  
    
    //loop over the rows
    for(int k=0;k<n-1;k++)
      {
	//find the pivot
	float_high_prec_t big=0.0;
	int ipiv=0;
	for(int i=k;i<n;i++) 
	  {
	    int ip=exch[i];
	    int ipk=n*ip+k;
	    float_high_prec_t size=abs(A[ipk])*x[ip];
	    if(size>big)
	      {
		big=size;
		ipiv=i;
	      }
	  }
	if(big.get_d()==0.0) crash("null big: %d",ipiv);
        
	std::swap(exch[ipiv],exch[k]);
	
	//pivotize
	float_high_prec_t pivot=A[n*exch[k]+k];
	for(int i=k+1;i<n;i++) 
	  {
	    A[n*exch[i]+k]/=pivot; //compute x to normalize other rows
	    for(int j=k+1;j<n;j++) //subtract from the other rows
	      A[n*exch[i]+j]-=A[n*exch[i]+k]*A[n*exch[k]+j];
	  }
      }
    if(A[n*exch[n-1]+n-1]==0.0) crash("last element null");
    
    //build solution
    for(int i=0;i<n;i++)
      {
	x[i]=b[exch[i]];
	for(int j=0;j<i;j++) x[i]-=A[n*exch[i]+j]*x[j];
      }
    
    //normalize with stored x
    for(int i=n-1;i>=0;i--) 
      {
	for(int j=i+1;j<n;j++) x[i]-=A[n*exch[i]+j]*x[j];
	x[i]/=A[n*exch[i]+i];
      }
  }
  
  //find nmax_err_points+1 extrema of Chebyshev polynomial
  void rat_approx_finder_t::find_cheb()
  {
    double coeff=(maximum-minimum)/(exp(1)-1);
    for(int i=0;i<nmax_err_points;i++)
      {
	xmax[i]=minimum+(exp(0.5*(1-cos((M_PI*i)/nzero_err_points)))-1)*coeff;
	zero[i]=minimum+(exp(0.5*(1-cos(M_PI*(2*i+1)/(2*nzero_err_points))))-1)*coeff;
      }
  }
  
  //set the steps
  void rat_approx_finder_t::set_step()
  {
    step[0]=zero[0]-minimum;
    for(int i=1;i<nzero_err_points;i++) step[i]=zero[i]-zero[i-1];
  }
  
  //perform an accomodation step
  void rat_approx_finder_t::new_step()
  {
    float_high_prec_t *yy=new float_high_prec_t[nmax_err_points];
    
    eclose=1e30;
    farther=0.0;
    
    //set left extrema
    float_high_prec_t zero0=minimum;
    
    for(int i=0;i<nmax_err_points;i++) 
      {
	//takes extrema
	float_high_prec_t zero1=zero[i],xm=xmax[i];
	if(i==nmax_err_points-1) zero1=maximum;
	
	//check if we can move in one of the dirs
	float_high_prec_t ym=get_abs_err(xm);
	float_high_prec_t q=step[i];
	float_high_prec_t xn=xm+q;
	float_high_prec_t yn;
	if(xn<zero0||!(xn<zero1)) // Cannot skip over adjacent boundaries
	  {
	    q=-q;
	    xn=xm;
	    yn=ym;
	  }
	else 
	  {
	    yn=get_abs_err(xn);
	    if(yn<ym)
	      {
		q=-q;
		xn=xm;
		yn=ym;
	      }
	  }
	
	//move until reaching barrier or decreasing error
	int istep=0,quit=0;
	while(!quit&&!(yn<ym))
	  {
	    istep++;
	    if(istep>9) quit=1;
	    else
	      {
		ym=yn;
		xm=xn;
		float_high_prec_t a=xm+q;
		if(a==xm||!(a>zero0)||!(a<zero1)) quit=1;
		else
		  {
		    xn=a;
		    yn=get_abs_err(xn);
		  }
	      }
	  }
	
	//copy new position and max
	xmax[i]=xm;
	yy[i]=ym;
	
	//search extream of error
	if(eclose>ym) eclose=ym;
	if(farther<ym) farther=ym;
	
	zero0=zero1;
      }
    
    verbosity_lv3_master_printf("eclose: %16.16lg, farther: %16.16lg, spread: %16.16lg, delta: %16.16lg\n",
		  eclose.get_d(),farther.get_d(),spread.get_d(),delta);
    
    //decrease step size if error spread increased
    float_high_prec_t q;
    //relative error spread
    if(eclose.get_d()!=0.0) q=farther/eclose-1;
    else q=farther;
    
    //spread is increasing: decrease step size
    if(!(q<spread)) delta*=0.5;
    
    spread=q;
    
    for(int i=0;i<nmax_err_points-1;i++) 
      {
	q=yy[i+1];
	if(q!=0.0) q=yy[i]/q-1;
	else q=0.0625;
	if(q>0.25) q=0.25;
	
	q*=xmax[i+1]-xmax[i];
	step[i]=q*delta;
      }
    step[nmax_err_points-1]=step[nmax_err_points-2];
    
    //insert new locations for the zeros
    for(int i=0;i<nzero_err_points;i++)
      {	
	float_high_prec_t xm=zero[i]-step[i];
	if(xm>minimum)
	  if(xm<maximum)
	    {
	      if(!(xm>xmax[i])) xm=0.5*(xmax[i]+zero[i]);
	      if(!(xm<xmax[i+1])) xm=0.5*(xmax[i+1]+zero[i]);
	      zero[i]=xm;
	    }
      }
    
    delete[] yy;
  }

  //set the linear system to be solved
  void rat_approx_finder_t::set_linear_system(float_high_prec_t *matr,float_high_prec_t *vec)
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

  //compute separately the numerator and denominator of the approximation
  void rat_approx_finder_t::compute_num_den_approx(float_high_prec_t &yn,float_high_prec_t &yd,float_high_prec_t x)
  {
    yn=coeff[degree];
    for(int i=degree-1;i>=0;i--) yn=x*yn+coeff[i];
    yd=x+coeff[2*degree];
    for(int i=2*degree-1;i>degree;i--) yd=x*yd+coeff[i];
  }
  
  //compute the approximation
  float_high_prec_t rat_approx_finder_t::compute_approx(float_high_prec_t x) 
  {
    //numerator and denominator
    float_high_prec_t yn,yd;
    compute_num_den_approx(yn,yd,x);
    return yn/yd;
  }
  
  //compute the error
  float_high_prec_t rat_approx_finder_t::get_err(float_high_prec_t x) 
  {
    //compute exact fun and approx
    float_high_prec_t fun=func_to_approx(x),app=compute_approx(x);
    
    //subtract
    float_high_prec_t err=fun-app;
    
    //normalize
    if(fun.get_d()!=0) err/=fun;
    
    return err;
  }
  
  //compute absolute relative error
  float_high_prec_t rat_approx_finder_t::get_abs_err(float_high_prec_t x)
  {return abs(get_err(x));}
  
  //calculate the roots of the approximation
  void rat_approx_finder_t::root_find(float_high_prec_t *roots,float_high_prec_t *poles,float_high_prec_t &cons) 
  {
    //define parameters
    const double upper=1,lower=-100000,tol=1.e-20;
    
    //find the numerator root
    float_high_prec_t *poly=new float_high_prec_t[nmax_err_points];
    for(int i=0;i<=degree;i++) poly[i]=coeff[i];
    
    for(int i=degree-1;i>=0;i--) 
      {
	roots[i]=root_find_Newton(poly,i+1,lower,upper,tol);
	if(roots[i]==0.0) crash("Failure to converge on root %ld/%d",i+1,degree);
	poly[0]/=-roots[i];
	for(int j=1;j<=i;j++) poly[j]=(poly[j-1]-poly[j])/roots[i];
      }
    
    //find the denominator roots
    poly[degree]=1;
    for(int i=0;i<degree;i++) poly[i]=coeff[degree+1+i];
    
    for(int i=degree-1;i>=0;i--) 
      {
	poles[i]=root_find_Newton(poly,i+1,lower,upper,tol);
	if(poles[i]==0.0) crash("Failure to converge on pole %ld/%d",i+1,degree);
	poly[0]/=-poles[i];
	for(int j=1;j<=i;j++) poly[j]=(poly[j-1]-poly[j])/poles[i];
      }
    
    cons=coeff[degree];
    
    delete[] poly;
  }
  
  //evaluate the polynomial
  float_high_prec_t rat_approx_finder_t::poly_eval(float_high_prec_t x,float_high_prec_t *poly,int size) 
  {
    float_high_prec_t f=poly[size];
    for(int i=size-1;i>=0;i--) f=f*x+poly[i];
    
    return f;
  }
  
  //evaluate the derivative
  float_high_prec_t rat_approx_finder_t::poly_der(float_high_prec_t x,float_high_prec_t *poly,int size) 
  {
    float_high_prec_t df=size*poly[size];
    for(int i=size-1;i>0;i--) df=df*x+i*poly[i];
    
    return df;
  }
  
  //Newton's method to calculate roots
  float_high_prec_t rat_approx_finder_t::root_find_Newton(float_high_prec_t *poly,int size,double x1,double x2,double acc) 
  {
    const int nmax_iter=10000;
    
    //start in the middle
    float_high_prec_t rtn=0.5*(x1+x2);
    
    //loop to find root
    int iter=0;
    float_high_prec_t dx;
    do
      {
	float_high_prec_t f=poly_eval(rtn,poly,size),df=poly_der(rtn,poly,size);
	dx=f/df;
	rtn-=dx;
	
	iter++;
      }
    while((iter<nmax_iter)&&(abs(dx)>=acc));
    
    if(iter==nmax_iter) crash("Maximum number of iterations exceeded in Newton_root");
    
    return rtn;
  }
  
  //decompose in a partial expansion
  void get_partial_fraction_expansion(float_high_prec_t *res,float_high_prec_t *poles,float_high_prec_t *roots,float_high_prec_t cons,int n)
  {
    float_high_prec_t *numerator=new float_high_prec_t[n],*denominator=new float_high_prec_t[n];
    for(int i=0;i<n;i++) res[i]=roots[i];
    
    //construct the polynomials explicitly
    numerator[0]=1.0;
    denominator[0]=1.0;
    for(int i=1;i<n;i++) numerator[i]=denominator[i]=0.0;
    
    for(int j=0;j<n;j++)
      for(int i=n-1;i>=0;i--)
	{
	  numerator[i]*=-res[j];
	  denominator[i]*=-poles[j];
	  
	  if(i>0)
	    {
	      numerator[i]+=numerator[i-1];
	      denominator[i]+=denominator[i-1];
	    }
	}
    
    //convert to proper fraction form, because now is in the form 1+n/d
    for(int i=0;i<n;i++) numerator[i]-=denominator[i];
    
    //find the residues of the partial fraction expansion and absorb the coefficients
    for(int i=0;i<n;i++)
      {
	res[i]=0.0;
	for(int j=n-1;j>=0;j--) res[i]=res[i]*poles[i]+numerator[j];
	
	for(int j=n-1;j>=0;j--) if(i!=j) res[i]=res[i]/(poles[i]-poles[j]);
	res[i]*=cons;
      }
    
    //res now holds the residues
    for(int i=0;i<n;i++) poles[i]=-poles[i];
      
    //move the ordering of the poles from smallest to largest
    for(int j=0;j<n;j++)
      {
	int small=j;
	for(int i=j+1;i<n;i++) if(poles[i]<poles[small]) small=i;
	
	if(small!=j)
	  {
	    std::swap(poles[small],poles[j]);
	    std::swap(res[small],res[j]);
	  }
      }
    
    delete[] numerator;
    delete[] denominator;
  }
  
  //generate the rational approximation
  double rat_approx_finder_t::generate_approx(float_high_prec_t *weights,float_high_prec_t *poles,float_high_prec_t &cons,double ext_minimum,double ext_maximum,int ext_degree,int ext_num,int ext_den,double target_err,double toll)
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
    step=new float_high_prec_t[nmax_err_points];
    coeff=new float_high_prec_t[nmax_err_points];
    zero=new float_high_prec_t[nmax_err_points];
    xmax=new float_high_prec_t[nmax_err_points];

    //set delta, initial spread and tolerance
    spread=1e37;
    delta=0.25;

    //#if HIGH_PREC==NATIVE_HIGH_PREC
    //approx_tolerance=2e-15;
    //#else
    //approx_tolerance=1e-15;
    //#endif
    approx_tolerance=toll;

    //set the initial guess and set step
    find_cheb();
    set_step();
    
    //iterate up to convergence
    int iter=0;
    float_high_prec_t *matr=new float_high_prec_t[nzero_err_points*nzero_err_points],*vec=new float_high_prec_t[nzero_err_points];
    do
      {
	// 1) set up the system to be solved
	set_linear_system(matr,vec);
	
	// 2) solve the system
	linear_system_solve(matr,coeff,vec,nzero_err_points);
	
	// 3) find maxima and minima
	if(iter==0 || (spread>approx_tolerance && farther>target_err)) new_step();
	if(delta<approx_tolerance)
	  {
	    master_printf("WARNING, reached precision %lg while computing %d terms approximation of x^(%d/%d) with tolerance %lg\n",
			  spread.get_d(),degree,num,den,approx_tolerance);
	    master_printf("precision not enough to reach %lg precision requested!!!\n",approx_tolerance);
#if HIGH_PREC==NATIVE_HIGH_PREC
	    master_printf("use GMP if possible!\n");
#else
	    master_printf("compile with higher precision!\n");
#endif
	  }
	iter++;
      }
    //while(float_high_prec_t_is_greater(spread,approx_tolerance));
    while(spread>approx_tolerance && delta>=approx_tolerance && eclose.get_d()<=target_err && farther.get_d()>target_err);
    
    //write some info
    if(spread<=approx_tolerance) verbosity_lv3_master_printf("Spread %lg reduced below %lg\n",spread.get_d(),approx_tolerance);
    if(eclose>target_err)  verbosity_lv3_master_printf("Accuracy cannot be better than %lg when %lg asked\n",eclose.get_d(),target_err);
    
    delete[] matr;
    delete[] vec;
    
    //get err at max and check
    if(farther.get_d()<=target_err) 
      {
	verbosity_lv2_master_printf("Converged in %d iters, maxerr %lg when asked %lg\n",iter,farther.get_d(),target_err);
	
	//compute the roots
	float_high_prec_t *roots=new float_high_prec_t[degree];
	root_find(roots,poles,cons);
	
	//decompose
	get_partial_fraction_expansion(weights,poles,roots,cons,degree);
	delete[] roots;
	
	for(int j=0;j<degree;j++)
	  verbosity_lv2_master_printf("Residue = %lg, Pole = %lg\n",weights[j].get_d(),poles[j].get_d());
	verbosity_lv2_master_printf("Const: %lg\n",cons.get_d());
      }
    else verbosity_lv2_master_printf("Not converged to %lg prec with %d poles in %d iters (reached: %lg)\n",target_err,degree,iter,farther.get_d());
      
    delete[] step;
    delete[] zero;
    delete[] coeff;
    delete[] xmax;
    
    //return the maximum error in the approximation
    return farther.get_d();
  }
  
  //generate an approximation
  double generate_approx(rat_approx_t &appr,double minimum,double maximum,int num,int den,double minerr,double tollerance)
  {
    appr.minimum=minimum;
    appr.maximum=maximum;
    appr.num=num;
    appr.den=den;
    
    //wrapper for 256bit output
    float_high_prec_t cons;
    float_high_prec_t *poles=new float_high_prec_t[appr.degree];
    float_high_prec_t *weights=new float_high_prec_t[appr.degree];
    
    //create the approx
    rat_approx_finder_t finder;
    double ans=finder.generate_approx(weights,poles,cons,minimum,maximum,appr.degree,num,den,minerr,tollerance);
    
    //copy
    appr.maxerr=ans;
    appr.cons=cons.get_d();
    for(int iterm=0;iterm<appr.degree;iterm++) appr.poles[iterm]=poles[iterm].get_d();
    for(int iterm=0;iterm<appr.degree;iterm++) appr.weights[iterm]=weights[iterm].get_d();
    
    delete[] poles;
    delete[] weights;
    
    return ans;
  }
  
  //generate an approximation
  void generate_approx_of_maxerr(rat_approx_t &appr,double minimum,double maximum,double maxerr,int num,int den,const char *name)
  {
    GET_THREAD_ID();
    
    //set the name of the approximation and deallocate poles
    if(name!=NULL) snprintf(appr.name,20,"%s",name);
    if(appr.degree!=0) rat_approx_destroy(&appr);
    
    //increase iteratively until it converges
    int degree=1;
    bool found=false;
    do
      {
	//allocate it
	rat_approx_create(&appr,degree,name);

	//generate
	double err;
	if(IS_MASTER_THREAD) err=generate_approx(appr,minimum,maximum,num,den,maxerr,0.01);
	THREAD_BROADCAST(err,err);
	
	//check if found
	found=(err<=maxerr);
	verbosity_lv3_master_printf("Approx x^(%x/%d) with %d poles can make an error of %lg when %lg required, found: %d\n",
	       num,den,degree,err,maxerr,found);
	
	//if not found increase number of poles
	THREAD_BARRIER();
	if(!found)
	  {
	    degree++;
	    rat_approx_destroy(&appr);
	  }
      }
    while(!found);
	
    //store required maximal error
    appr.maxerr=maxerr;    
  }
}
