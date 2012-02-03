#pragma once

//Compute the fermionic force the rooted staggered e/o improved theory.
//Passed conf do NOT have to contain the backfield.
//Of the result still need to be taken the TA
//The approximation need to be already scaled, and must contain physical mass term
void summ_the_sq4st_eoimpr_force(quad_su3 **F,quad_su3 **eo_conf,quad_u1 **u1b,color *pf,rat_approx *appr,double residue)
{
  //allocate each terms of the expansion
  //allocate temporary single solutions
  color *v_o[appr->nterms],*chi_e[appr->nterms];
  v_o[0]=nissa_malloc("v",(loc_volh+loc_bordh)*appr->nterms,color);
  chi_e[0]=nissa_malloc("chi",(loc_volh+loc_bordh)*appr->nterms,color);
  for(int iterm=1;iterm<appr->nterms;iterm++)
    {
      v_o[iterm]=v_o[iterm-1]+loc_volh+loc_bordh;
      chi_e[iterm]=chi_e[iterm-1]+loc_volh+loc_bordh;
    }
  
  //add the background fields
  add_backfield_to_conf(eo_conf,u1b);
  
  //invert the various terms
  inv_stD2ee_cgmm2s(chi_e,pf,eo_conf,appr->poles,appr->nterms,1000000,resiue,residue,0);
  
  //summ all the terms performing appropriate elaboration
  communicate_ev_color_borders(chi_e);
  for(int iterm=0;iterm<appr->nterm;iterm++)
    apply_st2Doe(v_0[iterm],eo+conf,chi_e[iterm]);
  
  //remove the background fields
  rem_backfield_from_conf(eo_conf,u1b);
  
  //conclude the calculation of the fermionic force
  memset(F,0,loc_vol*sizeof(quad_su3));
  for(int ivol=0;ivol<loc_volh;ivol++)
    {
      color w1[2],w2[2];
      for(int mu=0;mu<4;mu++)
	for(int iterm=0;iterm<appr->nterm;iterm++)
	  {
	    for(int ic1=0;ic1<3;ic1++)
	      for(int ic2=0;ic2<3;ic2++)
		{
		  complex temp;
		  
		  //this is for ivol=EVN
		  unsafe_complex_conj2_prod(temp,v_o[iterm][loceo_neighup[EVN][ivol][mu]][ic1],chi_e[iterm][ivol][ic2]);
		  //chpot to be added
		  summ_the_complex_prod_real(F[EVN][ivol][ic1][ic2],temp,appr->weight[iterm]);
		  
		  //this is for ivol=ODD
		  unsafe_complex_conj2_prod(temp,chi_e[iterm][loceo_neighup[ODD][ivol][mu]][ic1],v_o[iterm][ivol][ic2;//fin
		  //chpot to be added
		  subt_the_complex_prod_real(F[EVN][ivol][ic1][ic2],temp,appr->weight[iterm]);
	  }
    }
  
  //free the back field
  nissa_free(v_o);
  nissa_free(chi_e);
}
