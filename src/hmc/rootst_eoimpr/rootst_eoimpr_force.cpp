#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../geometry/geometry_eo.h"
#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../new_types/complex.h"
#include "../../inverters/staggered/cgm_invert_stD2ee_m2.h"
#include "../../linalgs/linalgs.h"
#include "../../operations/smearing/stout.h"

#include "../gauge/wilson_force.h"
#include "../gauge/tree_level_Symanzik_force.h"
#include "../backfield.h"

//Compute the fermionic force the rooted staggered e/o improved theory.
//Passed conf must NOT contain the backfield.
//Of the result still need to be taken the TA
//The approximation need to be already scaled, and must contain physical mass term
void summ_the_rootst_eoimpr_quark_force(quad_su3 **F,quad_su3 **eo_conf,color *pf,quad_u1 **u1b,rat_approx *appr,double residue)
{
  verbosity_lv1_master_printf("Computing quark force\n");
  
  //allocate each terms of the expansion
  color *v_o[appr->nterms],*chi_e[appr->nterms];
  for(int iterm=0;iterm<appr->nterms;iterm++)
    {
      v_o[iterm]=nissa_malloc("v_o",loc_volh+bord_volh,color);
      chi_e[iterm]=nissa_malloc("chi_e",loc_volh+bord_volh,color);
    }
  
  //add the background fields
  add_backfield_to_conf(eo_conf,u1b);
  
  //invert the various terms
  inv_stD2ee_m2_cgm_run_hm_up_to_mach_prec(chi_e,eo_conf,appr->poles,appr->nterms,1000000,residue,pf);
  
  //summ all the terms performing appropriate elaboration
  //possible improvement by communicating more borders together
  for(int iterm=0;iterm<appr->nterms;iterm++)
    apply_stDoe(v_o[iterm],eo_conf,chi_e[iterm]);
  
  //remove the background fields
  rem_backfield_from_conf(eo_conf,u1b);
  
  //communicate borders of v_o (could be improved...)
  for(int iterm=0;iterm<appr->nterms;iterm++) communicate_od_color_borders(v_o[iterm]);
  
  //conclude the calculation of the fermionic force
  nissa_loc_volh_loop(ivol)
    for(int mu=0;mu<4;mu++)
      for(int iterm=0;iterm<appr->nterms;iterm++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    {
	      //chpot to be added
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
  for(int iterm=0;iterm<appr->nterms;iterm++)
    {
      nissa_free(v_o[iterm]);
      nissa_free(chi_e[iterm]);
    }
}

//Finish the computation multiplying for the conf and taking TA
void full_rootst_eoimpr_force_finish_computation(quad_su3 **F,quad_su3 **conf)
{
  //remove the staggered phase from the conf, since it is already implemented in the force
  addrem_stagphases_to_eo_conf(conf);

  for(int eo=0;eo<2;eo++)
    nissa_loc_volh_loop(ivol)
      for(int mu=0;mu<4;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[eo][ivol][mu],F[eo][ivol][mu]);
	  unsafe_su3_traceless_anti_hermitian_part(F[eo][ivol][mu],temp);
	}
  
  //readd
  addrem_stagphases_to_eo_conf(conf);
}

//compute only the gauge part
void full_rootst_eoimpr_gluons_force(quad_su3 **F,quad_su3 **conf,theory_pars *physics)
{
  switch(physics->gac_type)
    {
    case Wilson_action: Wilson_force(F,conf,physics->beta);break;
    case tlSym_action: tree_level_Symanzik_force(F,conf,physics->beta);break;
    default: crash("Unknown action");
    }
  
  //add the stag phases to the force term, to cancel the one entering the force
  addrem_stagphases_to_eo_conf(F);

  full_rootst_eoimpr_force_finish_computation(F,conf);
}

//compute only the quark part, without stouting reampping
void full_rootst_eoimpr_quarks_force_no_stout_remapping(quad_su3 **F,quad_su3 **conf,color **pf,theory_pars *physics,rat_approx *appr,double residue)
{
  for(int eo=0;eo<2;eo++) vector_reset(F[eo]);
  for(int iflav=0;iflav<physics->nflavs;iflav++)
    summ_the_rootst_eoimpr_quark_force(F,conf,pf[iflav],physics->backfield[iflav],&(appr[iflav]),residue);
  
  //add the stag phases to the force term, coming from the disappered link in dS/d(U)
  addrem_stagphases_to_eo_conf(F);
}

//take into account the stout remapping procedure
void full_rootst_eoimpr_quarks_force(quad_su3 **F,quad_su3 **conf,color **pf,theory_pars *physics,rat_approx *appr,double residue)
{
  int nlev=physics->nstout_lev;
  
  //first of all we take care of the trivial case
  if(nlev==0) full_rootst_eoimpr_quarks_force_no_stout_remapping(F,conf,pf,physics,appr,residue);
  else
    {
      //allocate the stack of confs: conf is binded to sme_conf[0]
      quad_su3 ***sme_conf;
      stout_smear_conf_stack_allocate(sme_conf,conf,nlev);
      
      //smear iteratively retaining all the stack
      addrem_stagphases_to_eo_conf(sme_conf[0]); //remove the staggered phases
      stout_smear(sme_conf,conf,physics->stout_rho,nlev);
      
      //compute the force in terms of the most smeared conf
      addrem_stagphases_to_eo_conf(sme_conf[nlev]); //add to most smeared conf
      full_rootst_eoimpr_quarks_force_no_stout_remapping(F,sme_conf[nlev],pf,physics,appr,residue);
      
      //remap the force backward
      stouted_force_remap(F,sme_conf,physics->stout_rho,nlev);
      addrem_stagphases_to_eo_conf(sme_conf[0]); //add back again to the original conf
      
      //now free the stack of confs
      stout_smear_conf_stack_free(sme_conf,nlev);
    }
  
  full_rootst_eoimpr_force_finish_computation(F,conf);
}

//Compute the full force for the rooted staggered theory with nfl flavors
void full_rootst_eoimpr_force(quad_su3 **F,quad_su3 **conf,color **pf,theory_pars *physics,rat_approx *appr,double residue,hmc_force_piece piece)
{
  //First of all compute gluonic part of the force
  switch(piece)
    {
    case GAUGE_FORCE_ONLY: full_rootst_eoimpr_gluons_force(F,conf,physics);break;
    case QUARK_FORCE_ONLY: full_rootst_eoimpr_quarks_force(F,conf,pf,physics,appr,residue);break;
    case BOTH_FORCE_PIECES:
      //allocate an additional vector where to store fermionic pieces
      quad_su3 *F1[2]={nissa_malloc("F1[e]",loc_volh,quad_su3),nissa_malloc("F1[e]",loc_volh,quad_su3)};
      
      //compute the two pieces
      full_rootst_eoimpr_gluons_force(F,conf,physics);
      full_rootst_eoimpr_quarks_force(F1,conf,pf,physics,appr,residue);
  
      //summ the two pieces
      for(int eo=0;eo<2;eo++)
	{
	  double_vector_summassign((double*)F[eo],(double*)F1[eo],loc_volh*sizeof(quad_su3)/sizeof(double));
	  nissa_free(F1[eo]);
	}
      break;
    }
}
