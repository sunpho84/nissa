#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "base/debug.hpp"
#include "new_types/rat_approx.hpp"
#include "remez_algorithm.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  /*
  void linear_system_solve(double *A,double *x,double *b,int n)
  {
    double r[n],p[n],ap[n];
    for(int i=0;i<n;i++)
      {
	r[i]=p[i]=b[i];
	x[i]=0;
      }
    
    double rr=0;
    for(int i=0;i<n;i++) rr+=r[i]*r[i];
    double source_norm=0;
    for(int i=0;i<n;i++) source_norm+=b[i]*b[i];
    printf("Source_norm: %lg\n",source_norm);
    double rel_res;
    int iter=1;
    do
      {
	double pap=0;
	for(int i=0;i<n;i++)
	  {
	    ap[i]=0;
	    for(int j=0;j<n;j++) ap[i]+=A[i*n+j]*p[j];
	    pap+=ap[i]*p[i];
	  }
	//printf("pap: %lg\n",pap);
	
	double alpha=rr/pap;
	double roro=rr;
	
	//adjust new solution and residual,
	//compute new residual norm
	rr=0;
	for(int i=0;i<n;i++)
	  {
	    x[i]+=alpha*p[i];
	    r[i]-=alpha*ap[i];
	    rr+=r[i]*r[i]; //computing residual here we save one loop
	  }
	//printf("rr: %lg\n",rr);

	//adjust new krylov vector
	double beta=rr/roro;
	for(int i=0;i<n;i++) p[i]=r[i]+beta*p[i];
	
	//compute relative residue
	rel_res=(double)rr/source_norm;
	
	printf(" %d %lg\n",iter,rel_res);
	iter++;
      }
    while(rel_res>1e-29 && iter<100000);
  }

  void linear_system_solve_normal_equation(double *A,double *x,double *b,int n)
  {
    double A2[n*n];
    double bin[n];
    for(int i=0;i<n;i++)
      {
	bin[i]=0;
	for(int j=0;j<n;j++)
	  {
	    bin[i]+=A[j*n+i]*b[j];
	    A2[i*n+j]=0;
	    for(int k=0;k<n;k++) A2[i*n+j]+=A[k*n+i]*A[k*n+j];
	  }
      }
    linear_system_solve(A2,x,bin,n);
  }     
  */
  
  //solve the linear system A*x=b
  void linear_system_solve(std::vector<float_high_prec_t>& A,
			   std::vector<float_high_prec_t>& x,
			   const std::vector<float_high_prec_t>& b,
			   const int& n)
  {
    
    /*
    float_high_prec_t xcg[n];
    for(int i=0;i<n;i++) xcg[i]=0;
    
    double Ad[n*n],xd[n],bd[n];
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) Ad[i*n+j]=A[i*n+j].get_d();
    
    double source_norm=0;
    for(int i=0;i<n;i++) source_norm+=b[i].get_d()*b[i].get_d();
    
    double bdbd;
    int iter=0;
    double tol=sqr(pow(2.0,-high_prec_nbits()));
    do
      {
	bdbd=0;
	for(int i=0;i<n;i++)
	  {
	    float_high_prec_t t=b[i];
	    for(int j=0;j<n;j++) t-=A[i*n+j]*xcg[j];
	    bd[i]=t.get_d();
	    bdbd+=bd[i]*bd[i];
	  }
	bdbd/=source_norm;
	
	printf("accu res %d: %lg\n",iter,bdbd);
	
	if(bdbd>tol)
	  {
	    linear_system_solve_normal_equation(Ad,xd,bd,n);
	    for(int i=0;i<n;i++) xcg[i]+=xd[i];
	  }
	iter++;
      }
    while(bdbd>tol && iter<10000);
    
    for(int i=0;i<n;i++)
      printf("CG %d %lg\n",i,xcg[i].get_d());
      
    */
      
    std::vector<int> exch(n);
    
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
	if(row_norm==0.0) CRASH("num row norm");
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
	if(big.get_d()==0.0) CRASH("null big: %d",ipiv);
	
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
    if(A[n*exch[n-1]+n-1]==0.0) CRASH("last element null");
    
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
    
    /*
      for(int i=0;i<n;i++)
      printf("LU %d %lg\n",i,x[i].get_d());
      for(int i=0;i<n;i++)
      {
      float_high_prec_t d=x[i]-xcg[i];
      printf("DIFF %d %lg\n",i,d.get_d());
      }
      printf("Exp: %lg\n",pow(2.0,-high_prec_nbits()));
    */
    
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
    for(int i=1;i<nmax_err_points;i++)
      step[i]=zero[i]-zero[i-1];
  }
  
  //perform an accomodation step
  void rat_approx_finder_t::new_step(int iter)
  {
    
    std::vector<float_high_prec_t> yy(nmax_err_points);
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
    
    VERBOSITY_LV3_MASTER_PRINTF(" iter: %d, eclose: %16.16lg, farther: %16.16lg, spread: %16.16lg, delta: %16.16lg\n",
				iter,eclose.get_d(),farther.get_d(),spread.get_d(),delta);
    
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
  }
  
  //decompose in a partial expansion
  void get_partial_fraction_expansion(std::vector<float_high_prec_t>& res,
				      std::vector<float_high_prec_t>& poles,
				      const std::vector<float_high_prec_t>& roots,
				      const float_high_prec_t& cons,
				      const int& n)
  {
    std::vector<float_high_prec_t> numerator(n);
    std::vector<float_high_prec_t> denominator(n);
    
    for(int i=0;i<n;i++)
      res[i]=roots[i];
    
    //construct the polynomials explicitly
    numerator[0]=1.0;
    denominator[0]=1.0;
    for(int i=1;i<n;i++)
      numerator[i]=denominator[i]=0.0;
    
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
    for(int i=0;i<n;i++)
      poles[i]=-poles[i];
	
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
  }
  
  //generate an approximation
  double generate_approx(rat_approx_t& appr,
			 const double& minimum,
			 const double& maximum,
			 const int& num,
			 const int& den,
			 const double& minerr,
			 const double& tollerance)
  {
    
    appr.minimum=minimum;
    appr.maximum=maximum;
    appr.num=num;
    appr.den=den;
    
    //wrapper for 256bit output
    float_high_prec_t cons;
    std::vector<float_high_prec_t> poles(appr.degree());
    std::vector<float_high_prec_t> weights(appr.degree());
    
    //create the approx
    rat_approx_finder_t finder;
    const double ans=finder.generate_approx(weights,poles,cons,minimum,maximum,appr.degree(),num,den,minerr,tollerance);
    
    //copy
    appr.maxerr=ans;
    appr.cons=cons.get_d();
    for(int iterm=0;iterm<appr.degree();iterm++)
      appr.poles[iterm]=poles[iterm].get_d();
    for(int iterm=0;iterm<appr.degree();iterm++)
      appr.weights[iterm]=weights[iterm].get_d();
    
    return ans;
  }
  
  //generate an approximation
  void generate_approx_of_maxerr(rat_approx_t &appr,
				 const double& minimum,
				 const double& maximum,
				 const double& maxerr,
				 const int& num,
				 const int& den,
				 const char *name)
  {
    //perform a few checks
    if(num==den)
      CRASH("cannot work if the numerator has the same power of the denominator!");
    
    if(num==-den)
      {
	VERBOSITY_LV2_MASTER_PRINTF("Creating trivial approx for x^%d/%d\n",num,den);
	strncpy(appr.name,name,19);
	appr.resize(1);
	appr.num=num;
	appr.den=den;
	appr.maxerr=maxerr;
	appr.minimum=minimum;
	appr.maximum=maximum;
	appr.cons=0;
	appr.poles[0]=0;
	appr.weights[0]=1;
      }
    else
      {
        const double generate_time=take_time();
	
	//set the name of the approximation and deallocate poles
	if(name!=NULL) snprintf(appr.name,20,"%s",name);
	VERBOSITY_LV3_MASTER_PRINTF("Generating approximation of x^(%d/%d) with max error %lg over the interval [%lg,%lg]\n",
				    num,den,maxerr,minimum,maximum);
	
	//increase iteratively until it converges
	int degree=1;
	bool found=false;
	do
	  {
	    //allocate it
	    appr.resize(degree);
	    
	    //generate
	    double err;
	    err=generate_approx(appr,minimum,maximum,num,den,maxerr,0.01);
	    
	    //check if found
	    found=(err<=maxerr);
	    VERBOSITY_LV3_MASTER_PRINTF("Approx x^(%d/%d) with %d poles can make an error of %lg when %lg required, found: %d\n",
					num,den,degree,err,maxerr,found);
	    
	    //if not found increase number of poles
	    if(!found)degree++;
	  }
	while(!found);
	
	MASTER_PRINTF("Needed time: %lg s\n",take_time()-generate_time);
	
	//store required maximal error
	appr.maxerr=maxerr;
      }
  }
}
