#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "hmc/fermions/rootst_eoimpr_quark_force.hpp"
#include "hmc/fermions/roottm_clov_eoimpr_quark_force.hpp"
#include "hmc/hmc.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "hmc/theory_pars.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Finish the computation multiplying for the conf and taking TA
  void compute_quark_force_finish_computation(EoField<quad_su3>& F,
					      const EoField<quad_su3>& conf)
  {
    for(int eo=0;eo<2;eo++)
      {
	PAR(0,
	    locVolh,
	    CAPTURE(eo,
		    TO_WRITE(F),
		    TO_READ(conf)),
	    ivol,
	    {
	      for(int mu=0;mu<NDIM;mu++)
		{
		  su3 temp;
		  unsafe_su3_prod_su3(temp,conf[eo][ivol][mu],F[eo][ivol][mu]);
		  unsafe_su3_traceless_anti_hermitian_part(F[eo][ivol][mu],temp);
		}
	    });
      }
  }
  
  //compute the quark force, without stouting reampping
  void compute_quark_force_no_stout_remapping(EoField<quad_su3>& F,
					      EoField<quad_su3>& conf,
					      const std::vector<std::vector<pseudofermion_t>>& pf,
					      const theory_pars_t& tp,
					      const std::vector<rat_approx_t>& appr,
					      const double& residue)
  {
    //allocate or not clover term and inverse evn clover term
    EoField<clover_term_t>* Cl=nullptr;
    EvnField<inv_clover_term_t>* invCl_evn=nullptr;
    
    bool clover_to_be_computed=false;
    for(int iflav=0;iflav<tp.nflavs();iflav++)
      clover_to_be_computed|=ferm_discretiz::include_clover(tp.quarks[iflav].discretiz);
    if(clover_to_be_computed)
      {
	Cl=new EoField<clover_term_t>("Cl");
	invCl_evn=new EvnField<inv_clover_term_t>("inv_colver_term");
	chromo_operator(*Cl,conf);
      }
    
    //reset forces
    F.reset();
    
    for(int iflav=0;iflav<tp.nflavs();iflav++)
      {
	const quark_content_t &q=
	  tp.quarks[iflav];
	
	//if clover is included, compute it
	if(ferm_discretiz::include_clover(q.discretiz))
	  {
	    chromo_operator_include_cSW(*Cl,q.cSW);
	    invert_twisted_clover_term(*invCl_evn,q.mass,q.kappa,Cl->evenPart);
	  }
	
	const EoField<quad_u1>& bf=
	  tp.backfield[iflav];
	const rat_approx_t& app=
	  appr[iflav*nappr_per_quark+RAT_APPR_QUARK_FORCE];
	
	for(size_t ipf=0;ipf<pf[iflav].size();ipf++)
	  {
	    verbosity_lv2_master_printf("Computing quark force for flavour %d/%d, pseudofermion %zu/%zu\n",iflav+1,tp.nflavs(),ipf+1,pf[iflav].size());
	    
	    switch(q.discretiz)
	      {
	      case ferm_discretiz::ROOT_STAG:
		summ_the_rootst_eoimpr_quark_force(F,conf,*pf[iflav][ipf].stag,bf,app,residue);
		break;
	      case ferm_discretiz::ROOT_TM_CLOV:
		crash("reimplement");//summ_the_roottm_clov_eoimpr_quark_force(F,conf,q.kappa,q.cSW,Cl[ODD],invCl_evn,q.mass,*pf[iflav][ipf].Wils,bf,app,residue);
		break;
	      default:
		crash("unknown discretization");
	      }
	  }
	
	//remove cSW from chromo operator
	if(ferm_discretiz::include_clover(q.discretiz))
	  chromo_operator_remove_cSW(*Cl,q.cSW);
      }
    
    //free clover term if ever allocated
    if(clover_to_be_computed)
      {
	delete invCl_evn;
	delete Cl;
      }
  }
  
  //take into account the stout remapping procedure
  void compute_quark_force(EoField<quad_su3>& F,
			   EoField<quad_su3>& conf,
			   std::vector<std::vector<pseudofermion_t>>& pf,
			   const theory_pars_t& physics,
			   const std::vector<rat_approx_t>& appr,
			   const double& residue)
  {
    int nlevls=physics.stout_pars.nlevels;
    
    //first of all we take care of the trivial case
    if(nlevls==0)
      compute_quark_force_no_stout_remapping(F,conf,pf,physics,appr,residue);
    else
      {
	//allocate the stack of confs: conf is binded to sme_conf[0]
	std::vector<EoField<quad_su3>*> sme_conf;
	stout_smear_conf_stack_allocate(sme_conf,conf,nlevls);
	
	//smear iteratively retaining all the stack
	stout_smear_whole_stack(sme_conf,conf,physics.stout_pars);
	
	//compute the force in terms of the most smeared conf
	compute_quark_force_no_stout_remapping(F,*(sme_conf[nlevls]),pf,physics,appr,residue);
	
	//remap the force backward
	stouted_force_remap(F,sme_conf,physics.stout_pars);
	
	//now free the stack of confs
	stout_smear_conf_stack_free(sme_conf,nlevls);
      }
    
    compute_quark_force_finish_computation(F,conf);
    
    //print the intensity of the force
    // if(VERBOSITY_LV2)
    //   {
    // 	double norm=0;
    // 	for(int par=0;par<2;par++) norm+=double_vector_glb_norm2(F[par],locVolh);
    // 	master_printf("  Quark force average norm: %lg\n",sqrt(norm/glbVol));
    //   }
  }
}
