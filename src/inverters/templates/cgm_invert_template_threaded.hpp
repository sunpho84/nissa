#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/field.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  template <typename T,
	    typename F>
  void cgm_invert(std::vector<T>& sol,
		 const std::vector<double>& shifts,
		 F&& f,
		 const int& niter_max,
		 const std::vector<double>& residue,
		 const T& source)
  {
    const size_t nShift=
      shifts.size();
    
    const int n=
      source.nSites()*source.nInternalDegs;
    
    const int each=
      VERBOSITY_LV3?1:10;
    
    T s("s");
    T r("r");
    T p("p",WITH_HALO);
    std::vector<T> ps(nShift,"ps");
    
    //     -sol[*]=0
    //     -ps[*]=source
    for(size_t iShift=0;iShift<nShift;iShift++)
      {
	ps[iShift]=source;
	sol[iShift].reset();
      }
    
    //     -p=source
    //     -r=source
    //     -calculate source_norm=(r,r)
    p=source;
    r=source;
    const double source_norm2=r.norm2();
    
    //writes source norm
    verbosity_lv2_master_printf(" Source norm: %lg\n",source_norm2);
    if(source_norm2==0 || std::isnan(source_norm2)) crash("invalid norm: %lg",source_norm2);
    
    //writes initial residue
    verbosity_lv2_master_printf(" cgm iter 0 rel. residues: ");
    for(size_t iShift=0;iShift<nShift;iShift++)
      verbosity_lv2_master_printf("%1.4e  ",1.0);
    verbosity_lv2_master_printf("\n");
    
    int final_iter;
    double final_res[nShift];
    
    int iter=0;
    
    //     -betaa=1
    double betaa=1;
    
    //     -zps=zas=1
    //     -alphas=0
    double zps[nShift],zas[nShift],alphas[nShift];
    double zfs[nShift],betas[nShift];
    int run_flag[nShift],nrun_shift=nShift;
    for(int ishift=0;ishift<nShift;ishift++)
      {
	zps[ishift]=zas[ishift]=1;
	alphas[ishift]=0;
	run_flag[ishift]=1;
      }
    
    //     -alpha=0
    double alpha=0;
    
    //     -rr=(r,r)=source_norm
    double rr=source_norm2;
#ifdef CGM_DEBUG
    verbosity_lv3_master_printf("rr: %16.16lg\n",rr);
#endif
    
    double rfrf,betap;
    double res[nShift];
    
    do
      {
	//     this is already iteration 0
	final_iter=(++iter);
	
	//     -s=Ap
	if(use_async_communications and iter>1)
	  crash("reimplement");//CGM_FINISH_COMMUNICATING_BORDERS(p);
	
	cgm_inv_over_time+=take_time();
	f(s,0.0,p);
	cgm_inv_over_time-=take_time();
	
	//     -pap=(p,s)=(p,Ap)
	const double pap=
	  p.realPartOfScalarProdWith(s);
	
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("pap: %16.16lg (ap[0]: %16.16lg)\n",pap,((double*)s)[0]);
#endif
	//     calculate betaa=rr/pap=(r,r)/(p,Ap)
	betap=betaa;
	betaa=-rr/pap;
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("betap: %16.16lg, betaa: %16.16lg\n",betap,betaa);
#endif
	
	//     calculate
	//     -zfs
	//     -betas
	//     -x
	for(size_t ishift=0;ishift<nShift;ishift++)
	  {
	    if(run_flag[ishift]==1)
	      {
		const double ratio=
		  betap/(betaa*alpha*(1-zas[ishift]/zps[ishift])+betap*(1-shifts[ishift]*betaa));
		
		if(std::isnan(ratio))
		  crash("nanned");
		
		zfs[ishift]=zas[ishift]*ratio;
		betas[ishift]=betaa*ratio;
		
#ifdef CGM_DEBUG
		verbosity_lv3_master_printf("ishift %d [%lg] zas: %16.16lg, zps: %16.16lg, "
					    "zfs: %16.16lg, betas: %16.16lg\n",
					    ishift,shift[ishift],zas[ishift],zps[ishift],zfs[ishift],betas[ishift]);
#endif
		double_vector_summ_double_vector_prod_double((double*)(sol[ishift].data),(double*)(sol[ishift].data),(double*)(ps[ishift].data),-betas[ishift],n,DO_NOT_SET_FLAGS);
	      }
	    THREAD_BARRIER();
	  }
	
	//     calculate
	//     -r'=r+betaa*s=r+beta*Ap
	//     -rfrf=(r',r')
	double_vector_summ_double_vector_prod_double((double*)r.data,(double*)r.data,(double*)s.data,betaa,n);
	double_vector_glb_scalar_prod(&rfrf,(double*)r.data,(double*)r.data,n);
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("rfrf: %16.16lg\n",rfrf);
#endif
	
	//     calculate alpha=rfrf/rr=(r',r')/(r,r)
	alpha=rfrf/rr;
	
	//     calculate p'=r'+p*alpha
	double_vector_summ_double_vector_prod_double((double*)p.data,(double*)r.data,(double*)p.data,alpha,n);
	
	//start the communications of the border
	if(use_async_communications)
	  crash("reimplement");//CGM_START_COMMUNICATING_BORDERS(p);
	
	//     calculate
	//     -alphas=alpha*zfs*betas/zas*beta
	//     -ps'=r'+alpha*ps
	for(size_t ishift=0;ishift<nShift;ishift++)
	  {
	    if(run_flag[ishift]==1)
	      {
		alphas[ishift]=alpha*zfs[ishift]*betas[ishift]/(zas[ishift]*betaa);
#ifdef CGM_DEBUG
		verbosity_lv3_master_printf("ishift %d alpha: %16.16lg\n",ishift,alphas[ishift]);
#endif
		double_vector_linear_comb((double*)(ps[ishift].data),(double*)r.data,zfs[ishift],(double*)(ps[ishift].data),alphas[ishift],n,DO_NOT_SET_FLAGS);
		
		// shift z
		zps[ishift]=zas[ishift];
		zas[ishift]=zfs[ishift];
	      }
	    THREAD_BARRIER();
	  }
	
	//shift rr
	rr=rfrf;
	
	//check over residual
	if(iter%each==0)
	  verbosity_lv2_master_printf(" cgm iter %d rel. residues: ",iter);
	
	for(size_t ishift=0;ishift<nShift;ishift++)
	  if(run_flag[ishift])
	    {
	      final_res[ishift]=res[ishift]=rr*zfs[ishift]*zfs[ishift]/source_norm2;
	      if(iter%each==0)
		{
		  verbosity_lv2_master_printf("%1.4e  ",res[ishift]);
		  verbosity_lv3_master_printf("%16.16lg  ",res[ishift]);
		}
	      
	      if(res[ishift]<residue[ishift])
		{
		  run_flag[ishift]=0;
		  nrun_shift--;
		}
	    }
	  else
	    if(iter%each==0)
	      verbosity_lv2_master_printf(" * ");
	
	if(iter%each==0)
	  verbosity_lv2_master_printf("\n");
      }
    while(nrun_shift>0 and iter<niter_max);
    
    if(use_async_communications)
      crash("reimplement");//CGM_FINISH_COMMUNICATING_BORDERS(p);
    
    //print the final true residue
    for(int iShift=0;iShift<nShift;iShift++)
      {
	f(s,shifts[iShift],sol[iShift]);
	
	s-=source;
	const double res=s.norm2();
	
	verbosity_lv2_master_printf(" ishift %d, rel residue true=%lg approx=%lg commanded=%lg\n",
				    iShift,res/source_norm2,final_res[iShift],residue[iShift]);
	if(res/source_norm2>=2*final_res[iShift])
	  master_printf("WARNING: shift[%d]=%lg true residue (%lg) much larger than expected one (%lg)\n",
			iShift,shifts[iShift],res/source_norm2,final_res[iShift]);
      }
    
    verbosity_lv1_master_printf(" Total cgm iterations: %d\n",final_iter);
    
    //check if not converged
    if(final_iter==niter_max)
      crash("exit without converging");
    
    
    if(IS_MASTER_THREAD) cgm_inv_over_time+=take_time();
  }
}
