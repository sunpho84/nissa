#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <optional>

#include "base/bench.hpp"
#include "base/field.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  template <typename T,
	    typename F>
  void cg_invert(T& sol,
		 std::optional<T> guess,
		 F&& f,
		 const int& niter,
		 const double& residue,
		 const T& source)
  {
    verbosity_lv2_master_printf("\n");
    
    T s("s");
    T p("p",WITH_HALO);
    T r("r");
    
    //macro to be defined externally, allocating all the required additional vectors
    if(guess) sol=*guess;
    else sol.reset();
    
    START_TIMING(cg_inv_over_time,ncg_inv);
    int each=VERBOSITY_LV3?1:10;
    
    const double sourceNorm2=
      source.norm2();
    
    //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
    f(s,sol);
    
    r=source;
    r-=s;
    p=r;
    double delta=r.norm2();
    
    verbosity_lv2_master_printf("Source norm2: %lg\n",sourceNorm2);
    if(sourceNorm2==0 or std::isnan(sourceNorm2)) crash("invalid norm: %lg",sourceNorm2);
    verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/sourceNorm2);
    
    int final_iter;
    
    //main loop
    int iter=0;
    double lambda;
    do
      {
	//this is already iter 1
	final_iter=(++iter);
	
	//(r_k,r_k)/(p_k*DD*p_k)
	STOP_TIMING(cg_inv_over_time);
	f(s,p);
	cg_inv_over_time-=take_time();
	
	const double alpha=
	  s.realPartOfScalarProdWith(p);
	
	const double omega=
	  delta/alpha;
	
	//sol_(k+1)=x_k+omega*p_k
	sol.forEachSiteDeg([&p,omega](double& sol,const int& site,const int i)
	{
	  sol+=p(site,i)*omega;
	});
	
	//r_(k+1)=x_k-omega*p_k
	r.forEachSiteDeg([&s,omega](double& r,const int& site,const int i)
	{
	  r-=s(site,i)*omega;
	});
	
	//(r_(k+1),r_(k+1))
	lambda=r.norm2();
	
	//(r_(k+1),r_(k+1))/(r_k,r_k)
	const double gammag=
	  lambda/delta;
	delta=lambda;
	
	//checks
	if(std::isnan(gammag)) crash("nanned");
	
	//p_(k+1)=r_(k+1)+gammag*p_k
	p.forEachSiteDeg([&r,gammag](double& p,const int& site,const int i)
	{
	  p=r(site,i)+p*gammag;
	});
	
	if(iter%each==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,lambda/sourceNorm2);
      }
    while(lambda>=(residue*sourceNorm2) and iter<niter);
    
    //last calculation of residual
    f(s,sol);
    r=source;
    r-=s;
    lambda=r.norm2();
    
    verbosity_lv2_master_printf("final relative residue (after %d iters): %lg where %lg was required\n",
				final_iter,lambda/sourceNorm2,residue);
    if(lambda/sourceNorm2>=2*residue)
      master_printf("WARNING: true residue %lg much larger than required and expected one %lg\n",
		    lambda/sourceNorm2,residue);
    
    verbosity_lv1_master_printf(" Total cg iterations: %d\n",final_iter);
    
    //check if not converged
    if(final_iter==niter) crash("exit without converging");
    
    if(IS_MASTER_THREAD) cg_inv_over_time+=take_time();
  }
}
