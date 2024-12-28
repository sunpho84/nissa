#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/field.hpp"
#include "routines/ios.hpp"

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
    ncgm_inv++;
    cgm_inv_over_time-=take_time();
    
    const size_t nShift=
      shifts.size();
    
    const int each=
      VERBOSITY_LV3?1:10;
    
    T s("s");
    T t("t");
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
    const double source_norm2=
      r.norm2();
    
    //writes source norm
    verbosity_lv2_master_printf(" Source norm: %lg\n",source_norm2);
    if(source_norm2==0 or std::isnan(source_norm2)) crash("invalid norm: %lg",source_norm2);
    
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
    for(size_t ishift=0;ishift<nShift;ishift++)
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
    
	Field<double,T::fieldCoverage> buf("buf");
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
	cgm_inv_over_time1-=take_time();
	const double pap=
	  p.realPartOfScalarProdWith(s);
	cgm_inv_over_time1+=take_time();
	
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("pap: %16.16lg, ap[0]: %16.16lg, p[0]: %16.16lg\n",pap,s._data[0],p._data[0]);
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
		verbosity_lv3_master_printf("ishift %ld [%lg] zas: %16.16lg, zps: %16.16lg, "
					    "zfs: %16.16lg, betas: %16.16lg\n",
					    ishift,shifts[ishift],zas[ishift],zps[ishift],zfs[ishift],betas[ishift]);
#endif
		cgm_inv_over_time2-=take_time();
		FOR_EACH_SITE_DEG_OF_FIELD(r,
					   CAPTURE(sol=sol[ishift].getWritable(),
						   psi=ps[ishift].getReadable(),
						   beta=betas[ishift],
						   zfs=zfs[ishift]),
					   i,iD,
					   {
					     sol(i,iD)-=beta*psi(i,iD);
					   });
		cgm_inv_over_time2+=take_time();
	      }
	  }
	
	//     calculate
	//     -r'=r+betaa*s=r+beta*Ap
	//     -rfrf=(r',r')
	cgm_inv_over_time3-=take_time();
	FOR_EACH_SITE_DEG_OF_FIELD(r,
				   CAPTURE(TO_READ(s),
					   betaa,
					   TO_WRITE(r)),
				   i,iD,
				   {
				     r(i,iD)+=s(i,iD)*betaa;
				   });
	cgm_inv_over_time3+=take_time();
	
	cgm_inv_over_time4-=take_time();
	PAR(0,r.nSites(),
	  CAPTURE(TO_READ(r),
		  TO_WRITE(buf)),
	  site,
	  {
	    double s2{};
	    UNROLL_FOR(internalDeg,0,T::nInternalDegs)
	      s2+=sqr(t(site,internalDeg));
	    buf[site]=s2;
	  });
	cgm_inv_over_time4+=take_time();
	
#ifdef USE_CUDA
	cgm_inv_over_time5-=take_time();
	double rest=
	  thrust::reduce(thrust::device,
			 &buf[0],
			 &buf[r.nSites()],
			 0.0,
			 thrust::plus<double>());
	non_loc_reduce(&rest);
	cgm_inv_over_time5+=take_time();
#endif
	
	cgm_inv_over_time6-=take_time();
	glb_reduce(&res,buf,r.nSites());
	cgm_inv_over_time6+=take_time();
	
	cgm_inv_over_time7-=take_time();
	rfrf=r.norm2();
	cgm_inv_over_time7-=take_time();
	
#ifdef CGM_DEBUG
	verbosity_lv3_master_printf("rfrf: %16.16lg\n",rfrf);
#endif
	
	//     calculate alpha=rfrf/rr=(r',r')/(r,r)
	alpha=rfrf/rr;
	
	//     calculate p'=r'+p*alpha
	cgm_inv_over_time8-=take_time();
	FOR_EACH_SITE_DEG_OF_FIELD(r,
				   CAPTURE(TO_READ(r),
					   alpha,
					   TO_WRITE(p)),
				   i,iD,
				   {
				     p(i,iD)=r(i,iD)+p(i,iD)*alpha;
				   });
	cgm_inv_over_time8+=take_time();
	
	//start the communications of the border
	if(use_async_communications)
	  crash("reimplement");//CGM_START_COMMUNICATING_BORDERS(p);
	
	//     calculate
	//     -alphas=alpha*zfs*betas/zas*beta
	//     -ps'=zfs*r'+alpha*ps
	for(size_t ishift=0;ishift<nShift;ishift++)
	  {
	    if(run_flag[ishift]==1)
	      {
		alphas[ishift]=alpha*zfs[ishift]*betas[ishift]/(zas[ishift]*betaa);
#ifdef CGM_DEBUG
		verbosity_lv3_master_printf("ishift %ld alpha: %16.16lg\n",ishift,alphas[ishift]);
#endif
		cgm_inv_over_time9-=take_time();
		FOR_EACH_SITE_DEG_OF_FIELD(r,
					   CAPTURE(psi=ps[ishift].getWritable(),
						   TO_READ(r),
						   alpha=alphas[ishift],
						   zfs=zfs[ishift]),
					   i,iD,
					   {
					     psi(i,iD)=r(i,iD)*zfs+psi(i,iD)*alpha;
					   });
		cgm_inv_over_time9+=take_time();
		
		// shift z
		zps[ishift]=zas[ishift];
		zas[ishift]=zfs[ishift];
	      }
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
    for(size_t iShift=0;iShift<nShift;iShift++)
      {
	f(s,shifts[iShift],sol[iShift]);
	
	s-=source;
	const double res=
	  s.norm2();
	
	verbosity_lv2_master_printf(" ishift %zu, rel residue true=%lg approx=%lg commanded=%lg\n",
				    iShift,res/source_norm2,final_res[iShift],residue[iShift]);
	if(res/source_norm2>=2*final_res[iShift])
	  master_printf("WARNING: shift[%zu]=%lg true residue (%lg) much larger than expected one (%lg)\n",
			iShift,shifts[iShift],res/source_norm2,final_res[iShift]);
      }
    
    verbosity_lv1_master_printf(" Total cgm iterations: %d\n",final_iter);
    
    //check if not converged
    if(final_iter==niter_max)
      crash("exit without converging");
    
    cgm_inv_over_time+=take_time();
    
    master_printf("overhead 5 %lg s\n",cgm_inv_over_time5/ncgm_inv);
    master_printf("overhead 6 %lg s\n",cgm_inv_over_time6/ncgm_inv);
    master_printf("overhead 7 %lg s\n",cgm_inv_over_time7/ncgm_inv);
    master_printf("overhead 8 %lg s\n",cgm_inv_over_time8/ncgm_inv);
    master_printf("overhead 9 %lg s\n",cgm_inv_over_time9/ncgm_inv);
  }
}
