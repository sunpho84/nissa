#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "hmc/backfield.hpp"

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "inverters/staggered/cgm_32_invert_stD2ee_m2.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Compute the fermionic force the rooted staggered e/o improved theory.
  //Passed conf must NOT contain the backfield.
  //Of the result still need to be taken the TA
  //The approximation need to be already scaled, and must contain physical mass term
  THREADABLE_FUNCTION_8ARG(summ_the_rootst_eoimpr_quark_force, quad_su3**,F, double,charge, quad_su3**,eo_conf, color*,pf, int,quantization, quad_u1**,u1b, rat_approx_t*,appr, double,residue)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing quark force\n");
    
    //allocate each terms of the expansion
    color *v_o[appr->degree],*chi_e[appr->degree];
    for(int iterm=0;iterm<appr->degree;iterm++)
      {
	v_o[iterm]=nissa_malloc("v_o",loc_volh+bord_volh,color);
	chi_e[iterm]=nissa_malloc("chi_e",loc_volh+bord_volh,color);
      }
    
    //add the background fields
    add_backfield_to_conf(eo_conf,u1b);
    
    //invert the various terms
    inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(chi_e,eo_conf,appr->poles,appr->degree,1000000,residue,pf);
    
    ////////////////////
    
    
    //summ all the terms performing appropriate elaboration
    //possible improvement by communicating more borders together
    for(int iterm=0;iterm<appr->degree;iterm++) apply_stDoe(v_o[iterm],eo_conf,chi_e[iterm]);
        
    //remove the background fields
    rem_backfield_from_conf(eo_conf,u1b);
    
    //communicate borders of v_o (could be improved...)
    for(int iterm=0;iterm<appr->degree;iterm++) communicate_od_color_borders(v_o[iterm]);
    
    //conclude the calculation of the fermionic force
    for(int iterm=0;iterm<appr->degree;iterm++)
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	for(int mu=0;mu<NDIM;mu++)
	  for(int ic1=0;ic1<NCOL;ic1++)
	    for(int ic2=0;ic2<NCOL;ic2++)
	      {
		complex temp1,temp2;
		
		//this is for ivol=EVN
		unsafe_complex_conj2_prod(temp1,v_o[iterm][loceo_neighup[EVN][ivol][mu]][ic1],chi_e[iterm][ivol][ic2]);
		unsafe_complex_prod(temp2,temp1,u1b[EVN][ivol][mu]);
		complex_summ_the_prod_double(F[EVN][ivol][mu][ic1][ic2],temp2,appr->weights[iterm]);
		
		//this is for ivol=ODD
		unsafe_complex_conj2_prod(temp1,chi_e[iterm][loceo_neighup[ODD][ivol][mu]][ic1],v_o[iterm][ivol][ic2]);
		unsafe_complex_prod(temp2,temp1,u1b[ODD][ivol][mu]);
		complex_subt_the_prod_double(F[ODD][ivol][mu][ic1][ic2],temp2,appr->weights[iterm]);
	      }
    
    //free
    for(int iterm=0;iterm<appr->degree;iterm++)
      {
	nissa_free(v_o[iterm]);
	nissa_free(chi_e[iterm]);
      }
  }
  THREADABLE_FUNCTION_END

  //Finish the computation multiplying for the conf and taking TA
  THREADABLE_FUNCTION_2ARG(compute_rootst_eoimpr_quark_force_finish_computation, quad_su3**,F, quad_su3**,conf)
  {
    GET_THREAD_ID();
    
    //remove the staggered phase from the conf, since they are already implemented in the force
    addrem_stagphases_to_eo_conf(conf);
    
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
    
    //readd
    addrem_stagphases_to_eo_conf(conf);
  }
  THREADABLE_FUNCTION_END
  
  //compute the quark force, without stouting reampping
  THREADABLE_FUNCTION_7ARG(compute_rootst_eoimpr_quark_force_no_stout_remapping, quad_su3**,F, quad_su3**,conf, color***,pf, theory_pars_t*,tp, rat_approx_t*,appr, int*,npfs, double,residue)
  {
    //reset forces
    for(int eo=0;eo<2;eo++) vector_reset(F[eo]);
    
    for(int iflav=0;iflav<tp->nflavs;iflav++)
      for(int ipf=0;ipf<npfs[iflav];ipf++)
	summ_the_rootst_eoimpr_quark_force(F,tp->quark_content[iflav].charge,conf,pf[iflav][ipf],tp->em_field_pars.flag,tp->backfield[iflav],appr+(iflav*3+2),residue); //flag set quantization
    
    //add the stag phases to the force term, coming from the disappered link in dS/d(U)
    addrem_stagphases_to_eo_conf(F);
  }
  THREADABLE_FUNCTION_END
  
  //take into account the stout remapping procedure
  THREADABLE_FUNCTION_7ARG(compute_rootst_eoimpr_quark_force, quad_su3**,F, quad_su3**,conf, color***,pf, theory_pars_t*,physics, rat_approx_t*,appr, int*,npfs, double,residue)
  {
    int nlevls=physics->stout_pars.nlevels;
    
    //first of all we take care of the trivial case
    if(nlevls==0) compute_rootst_eoimpr_quark_force_no_stout_remapping(F,conf,pf,physics,appr,npfs,residue);
    else
      {
	//allocate the stack of confs: conf is binded to sme_conf[0]
	quad_su3 ***sme_conf;
	stout_smear_conf_stack_allocate(&sme_conf,conf,nlevls);
	
	//smear iteratively retaining all the stack
	addrem_stagphases_to_eo_conf(sme_conf[0]); //remove the staggered phases
	stout_smear_whole_stack(sme_conf,conf,&(physics->stout_pars));
	
	//compute the force in terms of the most smeared conf
	addrem_stagphases_to_eo_conf(sme_conf[nlevls]); //add to most smeared conf
	compute_rootst_eoimpr_quark_force_no_stout_remapping(F,sme_conf[nlevls],pf,physics,appr,npfs,residue);
	
	//remap the force backward
	stouted_force_remap(F,sme_conf,&(physics->stout_pars));
	addrem_stagphases_to_eo_conf(sme_conf[0]); //add back again to the original conf
	
	//now free the stack of confs
	stout_smear_conf_stack_free(&sme_conf,nlevls);
      }
    
    compute_rootst_eoimpr_quark_force_finish_computation(F,conf);
  }
  THREADABLE_FUNCTION_END
}
