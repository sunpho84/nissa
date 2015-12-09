#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "hmc/rootst_eoimpr/rootst_eoimpr_quark_force.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Finish the computation multiplying for the conf and taking TA
  THREADABLE_FUNCTION_2ARG(compute_quark_force_finish_computation, quad_su3**,F, quad_su3**,conf)
  {
    GET_THREAD_ID();
    
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3 temp;
	      unsafe_su3_prod_su3(temp,conf[eo][ivol][mu],F[eo][ivol][mu]);
	      unsafe_su3_traceless_anti_hermitian_part(F[eo][ivol][mu],temp);
	    }
      }
  }
  THREADABLE_FUNCTION_END
  
  //compute the quark force, without stouting reampping
  THREADABLE_FUNCTION_7ARG(compute_quark_force_no_stout_remapping, quad_su3**,F, quad_su3**,conf, pseudofermion_t*,pf, theory_pars_t*,tp, rat_approx_t*,appr, int*,npfs, double,residue)
  {
    //reset forces
    for(int eo=0;eo<2;eo++) vector_reset(F[eo]);
    
    for(int iflav=0;iflav<tp->nflavs;iflav++)
      for(int ipf=0;ipf<npfs[iflav];ipf++)
	if(tp->quark_content[iflav].is_stag)
	summ_the_rootst_eoimpr_quark_force(F,tp->quark_content[iflav].charge,conf,pf[iflav].stag[ipf],tp->em_field_pars.flag,tp->backfield[iflav],appr+(iflav*3+2),residue);
	else crash("non staggered not yet implemented");
    
  }
  THREADABLE_FUNCTION_END
  
  //take into account the stout remapping procedure
  THREADABLE_FUNCTION_7ARG(compute_quark_force, quad_su3**,F, quad_su3**,conf, pseudofermion_t*,pf, theory_pars_t*,physics, rat_approx_t*,appr, int*,npfs, double,residue)
  {
    int nlevls=physics->stout_pars.nlevels;
    
    //first of all we take care of the trivial case
    if(nlevls==0) compute_quark_force_no_stout_remapping(F,conf,pf,physics,appr,npfs,residue);
    else
      {
	//allocate the stack of confs: conf is binded to sme_conf[0]
	quad_su3 ***sme_conf;
	stout_smear_conf_stack_allocate(&sme_conf,conf,nlevls);
	
	//smear iteratively retaining all the stack
	stout_smear_whole_stack(sme_conf,conf,&(physics->stout_pars));
	
	//compute the force in terms of the most smeared conf
	compute_quark_force_no_stout_remapping(F,sme_conf[nlevls],pf,physics,appr,npfs,residue);
	
	//remap the force backward
	stouted_force_remap(F,sme_conf,&(physics->stout_pars));
	
	//now free the stack of confs
	stout_smear_conf_stack_free(&sme_conf,nlevls);
      }
    
    compute_quark_force_finish_computation(F,conf);
  }
  THREADABLE_FUNCTION_END
}
