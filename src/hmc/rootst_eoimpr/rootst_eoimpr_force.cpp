#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../new_types/complex.h"
#include "../../inverters/staggered/cgm_invert_stD2ee_m2.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"

#include "../gauge/wilson_force.h"
#include "../backfield.h"

//Compute the fermionic force the rooted staggered e/o improved theory.
//Passed conf must NOT contain the backfield.
//Of the result still need to be taken the TA
//The approximation need to be already scaled, and must contain physical mass term
void summ_the_rootst_eoimpr_quarks_force(quad_su3 **F,quad_su3 **eo_conf,color *pf,quad_u1 **u1b,rat_approx *appr,double residue)
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

//Compute the full force for the rooted staggered theory with nfl flavors
//Conf and pf borders are assumed to have been already communicated
void full_rootst_eoimpr_force(quad_su3 **F,quad_su3 **conf,color **pf,theory_pars *physic,rat_approx *appr,double residue,hmc_force_piece piece)
{
  //First of all compute gluonic part of the force
  if(piece==GAUGE_ONLY_FORCE||piece==BOTH_FORCE_PARTS)
    wilson_force(F,conf,physic->beta);
  else
    for(int eo=0;eo<2;eo++)
      vector_reset(F[eo]);

  //Then summ the contributes coming from each flavor
  if(piece==QUARK_ONLY_FORCE||piece==BOTH_FORCE_PARTS)
    for(int iflav=0;iflav<physic->nflavs;iflav++)
      summ_the_rootst_eoimpr_quarks_force(F,conf,pf[iflav],physic->backfield[iflav],&(appr[iflav]),residue);
  
  //Finish the computation multiplying for the conf and taking TA
  for(int eo=0;eo<2;eo++)
    nissa_loc_volh_loop(ivol)
      for(int mu=0;mu<4;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[eo][ivol][mu],F[eo][ivol][mu]);
	  su3_traceless_anti_hermitian_part(F[eo][ivol][mu],temp);
	}
}
