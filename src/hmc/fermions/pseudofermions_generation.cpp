#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "hmc/hmc.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  // Generate the pseudofermions in the root tm case
  void generate_root_tm_clov_pseudo_fermion(EvnField<spincolor>& pf,
					    const EoField<quad_su3>& conf,
					    const EoField<quad_u1>& u1b,
					    const rat_approx_t& rat,
					    const double& residue,
					    const quark_content_t& q,
					    const EvnField<spincolor>& eta)
  {
    CRASH("reimplement");
    // //allocate and compute clover term
    // EoField<clover_term_t> Cl("Cl");
    // EoField<inv_clover_term_t> invCl_evn("invCl_evn");
    // chromo_operator(Cl,conf);
    
    // chromo_operator_include_cSW(Cl,q.cSW);
    // invert_twisted_clover_term(invCl_evn,q.mass,q.kappa,Cl[EVN]);
    
    // add_backfield_without_stagphases_to_conf(conf,u1b);
    // if(2*rat->num==rat->den)
    //   {
    // 	spincolor *tmp=nissa_malloc("tmp1",locVolh+bord_volh,spincolor);
    // 	tmclovDkern_eoprec_eos(pf,tmp,conf,q.kappa,Cl[ODD],invCl_evn,false,q.mass,eta);
    // 	nissa_free(tmp);
    //   }
    // else
    //   summ_src_and_all_inv_tmclovDkern_eoprec_square_portable(pf,conf,q.kappa,Cl[ODD],invCl_evn,q.mass,rat,10000000,residue,eta);
    // rem_backfield_without_stagphases_from_conf(conf,u1b);
    
    // //free clover term
    // nissa_free(invCl_evn);
    // for(int eo=0;eo<2;eo++) nissa_free(Cl[eo]);
  }
  
  // Generate the psewudoferimions in the root stag case
  void generate_root_stag_pseudo_fermion(EvnField<color>& pf,
					 EoField<quad_su3>& conf,
					 const EoField<quad_u1>& u1b,
					 const rat_approx_t& rat,
					 const double& residue,
					 const quark_content_t& q,
					 const EvnField<color>& eta)
  {
    add_backfield_with_stagphases_to_conf(conf,u1b);
    summ_src_and_all_inv_stD2ee_m2_cgm(pf,conf,rat,10000000,residue,eta);
    rem_backfield_with_stagphases_from_conf(conf,u1b);
  }
  
  //generate pseudo-fermion
  double generate_pseudo_fermion(pseudofermion_t& pf,
				 EoField<quad_su3>& conf,
				 const EoField<quad_u1>& u1b,
				 const rat_approx_t& rat,
				 const double& residue,
				 const quark_content_t& q)
  {
    //generate the random field
    pseudofermion_t pf_hb_vec(q.discretiz);
    pf_hb_vec.fill();
    
    //invert to perform hv
    switch(q.discretiz)
      {
      case ferm_discretiz::ROOT_STAG:
	generate_root_stag_pseudo_fermion(*pf.stag,conf,u1b,rat,residue,q,*pf_hb_vec.stag);
	break;
      case ferm_discretiz::ROOT_TM_CLOV:
	generate_root_tm_clov_pseudo_fermion(*pf.Wils,conf,u1b,rat,residue,q,*pf_hb_vec.Wils);
	break;
      default:
	CRASH("not supported");
	break;
      }
    
    //compute action
    const double action=
      pf_hb_vec.norm2();
    
    return action;
  }
  
  //gemerate all pseudofermions
  double generate_pseudofermions(std::vector<std::vector<pseudofermion_t>>& pf,
				 EoField<quad_su3>& conf,
				 const theory_pars_t& theory_pars,
				 const hmc_evol_pars_t& simul_pars,
				 const std::vector<rat_approx_t>& rat_appr)
  {
    //create pseudo-fermions and store action
    double pf_action=0;
    for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
      {
	//get quark
	const quark_content_t &q=theory_pars.quarks[iflav];
	
	//loop over pseudofermions
	for(int ipf=0;ipf<simul_pars.npseudo_fs[iflav];ipf++)
	  {
	    VERBOSITY_LV1_MASTER_PRINTF("Generating pseudofermion %d/%d for flavour %d/%d\n",ipf+1,simul_pars.npseudo_fs[iflav],iflav+1,theory_pars.nflavs());
	    const double pf_action_flav=
	      generate_pseudo_fermion(pf[iflav][ipf],conf,theory_pars.backfield[iflav],rat_appr[nappr_per_quark*iflav+RAT_APPR_PF_GEN],simul_pars.pf_action_residue,q);
	    pf_action+=pf_action_flav;
	  }
      }
    
    return pf_action;
  }
}
