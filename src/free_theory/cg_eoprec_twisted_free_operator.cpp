#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <optional>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"
#include "free_theory_types.hpp"
#include "twisted_free_Dirac_eoprec_operator.hpp"

namespace nissa
{
  //invert Koo defined in equation (7)
  void inv_tmDkern_eoprec_square_eos(OddField<spin>& sol,
				     const std::optional<OddField<spin>>& guess,
				     const tm_quark_info& qu,
				     const int& nMaxIter,
				     const double& residue,
				     const OddField<spin>& source)
  {
    crash("reimplement");
    
    // int niter=nMaxIter;
    // int riter=0;
    // int rniter=5;
    // OddField<spin> p("p",WITH_HALO);
    // OddField<spin> r("r"),s("s");
    // OddField<spin> temp1("temp1",WITH_HALO);
    // EvnField<spin> temp2("temp2",WITH_HALO);
    
    // ///////////////// prepare the internal source /////////////////
    
    // if(guess) sol=*guess;
    // else sol.reset();
    
    // //external loop, used if the internal exceed the maximal number of iterations
    // double lambda; //(r_(k+1),r_(k+1))
    // const double sourceNorm2=source.norm2();
    // do
    //   {
    // 	//calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
    // 	tmDkern_eoprec_square_eos(s,temp1,temp2,qu,sol);
	
    // 	FOR_EACH_SITE_DEG_OF_FIELD(r,
    // 				   CAPTURE(TO_WRITE(r),
    // 					   TO_READ(source),
    // 					   TO_READ(s)),site,i,
    // 	{
    // 	  r(site,i)=source(site,i)-s(site,i);
    // 	});
    // 	p=r;
    // 	double delta=r.norm2();
	
    // 	if(riter==0)
    // 	  {
    // 	    master_printf("\nSource norm: %lg\n",sourceNorm2);
    // 	    master_printf("iter 0 relative residue: %lg\n",delta/sourceNorm2);
    // 	  }
	
    // 	//main loop
    // 	int iter=0;
    // 	do
    // 	  {
    // 	    tmDkern_eoprec_square_eos(s,temp1,temp2,qu,p);
    // 	    const double alpha=s.realPartOfScalarProdWith(p);
    // 	    const double omega=delta/alpha;
	    
    // 	    FOR_EACH_SITE_DEG_OF_FIELD(sol,
    // 				       CAPTURE(TO_WRITE(sol),
    // 					       TO_READ(p),
    // 					       omega),site,i,
    // 	    {
    // 	      sol(site,i)+=p(site,i)*omega;
    // 	    });
    // 	    r.forEachSiteDeg([&s,omega](double& r,const int& site,const int i)
    // 	    {
    // 	      r-=s(site,i)*omega;
    // 	    });
    // 	    lambda=r.norm2();
	    
    // 	    const double gammag=lambda/delta;
    // 	    delta=lambda;
	    
    // 	    //p_(k+1)=r_(k+1)+gammag*p_k
    // 	    p.forEachSiteDeg([&r,gammag](double& p,const int& site,const int i)
    // 	    {
    // 	      p=r(site,i)+p*gammag;
    // 	    });
	    
    // 	    iter++;
	    
    // 	    if(iter%10==0)
    // 	      master_printf("iter %d relative residue: %lg\n",iter,lambda/sourceNorm2);
    // 	  }
    // 	while(lambda>(residue*sourceNorm2) and iter<niter);
	
    // 	//last calculation of residual, in the case iter>niter
    // 	tmDkern_eoprec_square_eos(s,temp1,temp2,qu,sol);
    // 	r.forEachSiteDeg([&s,&source](double& r,const int& site,const int i)
    // 	{
    // 	  r=source(site,i)-s(site,i);
    // 	});
    // 	lambda=r.norm2();
	
    // 	master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",iter,lambda/sourceNorm2,residue);
	
    // 	riter++;
    //   }
    // while(lambda>(residue*sourceNorm2) && riter<rniter);
  }
  
  //Invert twisted mass operator using e/o preconditioning.
  void inv_tmD_cg_eoprec_eos(LxField<spin>& solution_lx,
			     std::optional<OddField<spin>> guess_Koo,
			     const tm_quark_info& qu,
			     const int& nitermax,
			     const double& residue,
			     const LxField<spin>& source_lx)
  {
    crash("reimplement");
    
    // //prepare the e/o split version of the source
    // EoField<spin> source_eos("source_eos",WITH_HALO);
    // split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    // //prepare the e/o split version of the solution
    // EoField<spin> solution_eos("solution_eos",WITH_HALO);
    
    // ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
    
    // OddField<spin> varphi("varphi",WITH_HALO);
    
    // //Equation (8.a)
    // OddField<spin> temp("temp",WITH_HALO);
    // inv_tmDee_or_oo_eos(temp,qu,source_eos[EVN]);
    
    // //Equation (8.b)
    // tmn2Deo_or_tmn2Doe_eos(varphi,temp,qu.bc);
    // NISSA_PARALLEL_LOOP(ivol,0,locVolh)
    //   for(int id=0;id<2;id++)
    // 	for(int ri=0;ri<2;ri++)
    // 	  { //gamma5 is explicitely wrote
    // 	    varphi[ivol][id  ][ri]=+source_eos[ODD][ivol][id  ][ri]+varphi[ivol][id  ][ri]*0.5;
    // 	    varphi[ivol][id+2][ri]=-source_eos[ODD][ivol][id+2][ri]-varphi[ivol][id+2][ri]*0.5;
    // 	  }
    // NISSA_PARALLEL_LOOP_END;
    // varphi.invalidateHalo();
    
    // //Equation (9)
    // inv_tmDkern_eoprec_square_eos(temp,guess_Koo,qu,nitermax,residue,varphi);
    // tm_quark_info mqu=qu;
    // mqu.mass*=-1;
    // //auto& uuu=solution_eos[EVN].castSitesCoverage<>()
    // tmDkern_eoprec_eos(solution_eos.oddPart,solution_eos.evenPart,mqu,temp);
    // if(guess_Koo) (*guess_Koo)=temp; //if a guess was passed, return new one
    
    // //Equation (10)
    // tmn2Deo_or_tmn2Doe_eos(varphi,solution_eos.oddPart,qu.bc);
    // NISSA_PARALLEL_LOOP(ivol,0,locVolh)
    //   for(int id=0;id<NDIRAC;id++)
    // 	for(int ri=0;ri<2;ri++)
    // 	  varphi[ivol][id][ri]=source_eos[EVN][ivol][id][ri]+varphi[ivol][id][ri]*0.5;
    // NISSA_PARALLEL_LOOP_END;
    // varphi.invalidateHalo();
    // inv_tmDee_or_oo_eos(solution_eos.evenPart,qu,varphi);
    
    // /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    // paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
  }
}
