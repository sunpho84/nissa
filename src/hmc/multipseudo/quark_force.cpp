#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "hmc/fermions/rootst_eoimpr_quark_force.hpp"
#include "hmc/hmc.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "hmc/theory_pars.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Finish the computation multiplying for the conf and taking TA
  THREADABLE_FUNCTION_2ARG(compute_quark_force_finish_computation, eo_ptr<quad_su3>,F, eo_ptr<quad_su3>,conf)
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
	NISSA_PARALLEL_LOOP_END;
	
      }
  }
  THREADABLE_FUNCTION_END
  
  //compute the quark force, without stouting reampping
  THREADABLE_FUNCTION_6ARG(compute_quark_force_no_stout_remapping, quad_su3**,F, quad_su3**,conf, std::vector<std::vector<pseudofermion_t> >*,pf, theory_pars_t*,tp, std::vector<rat_approx_t>*,appr, double,residue)
  {
    //reset forces
    for(int eo=0;eo<2;eo++) vector_reset(F[eo]);
    
    for(int iflav=0;iflav<tp->nflavs();iflav++)
      for(size_t ipf=0;ipf<(*pf)[iflav].size();ipf++)
	{
	  verbosity_lv2_master_printf("Computing quark force for flavour %d/%d, pseudofermion %d/%d\n",iflav+1,tp->nflavs(),ipf+1,(*pf)[iflav].size());
	  
	  switch(tp->quarks[iflav].discretiz)
	    {
	    case ferm_discretiz::ROOT_STAG:
	      summ_the_rootst_eoimpr_quark_force(F,tp->quarks[iflav].charge,conf,(*pf)[iflav][ipf].stag,tp->em_field_pars.flag,tp->backfield[iflav],&((*appr)[iflav*nappr_per_quark+RAT_APPR_QUARK_FORCE]),residue);break;
	    default:
	      crash("non staggered not yet implemented");
	    }
	}
    
  }
  THREADABLE_FUNCTION_END
  
  //take into account the stout remapping procedure
  THREADABLE_FUNCTION_6ARG(compute_quark_force, quad_su3**,F, quad_su3**,conf, std::vector<std::vector<pseudofermion_t> >*,pf, theory_pars_t*,physics, std::vector<rat_approx_t>*,appr, double,residue)
  {
    int nlevls=physics->stout_pars.nlevels;
    
    //first of all we take care of the trivial case
    if(nlevls==0) compute_quark_force_no_stout_remapping(F,conf,pf,physics,appr,residue);
    else
      {
	//allocate the stack of confs: conf is binded to sme_conf[0]
	quad_su3 ***sme_conf;
	stout_smear_conf_stack_allocate(&sme_conf,conf,nlevls);
	
	//smear iteratively retaining all the stack
	stout_smear_whole_stack(sme_conf,conf,&(physics->stout_pars));
	
	//compute the force in terms of the most smeared conf
	compute_quark_force_no_stout_remapping(F,sme_conf[nlevls],pf,physics,appr,residue);
	
	//remap the force backward
	stouted_force_remap(F,sme_conf,&(physics->stout_pars));
	
	//now free the stack of confs
	stout_smear_conf_stack_free(&sme_conf,nlevls);
      }
    
    compute_quark_force_finish_computation(F,conf);
    
    //print the intensity of the force
    if(VERBOSITY_LV2)
      {
	double norm=0;
	for(int par=0;par<2;par++) norm+=double_vector_glb_norm2(F[par],loc_volh);
	master_printf("  Quark force average norm: %lg\n",sqrt(norm/glb_vol));
      }
  }
  THREADABLE_FUNCTION_END
}
