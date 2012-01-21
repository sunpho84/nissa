#pragma once

//Compute the fermionic force the rooted staggered e/o improved theory.
//Passed conf do NOT have to contain the backfield.
//Of the result still need to be taken the TA
//The approximation need to be already scaled, and must contain physical mass term
void summ_the_r4st_eoimpr_force(quad_su3 *F,quad_su3 **eo_conf,quad_u1 **u1b,color *pf,rat_approx *appr,double residue)
{
  //add the background fields
  add_backfield_to_conf(eo_conf,u1b);
  
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
  
  //invert the various terms
  inv_stD2ee_cgmm2s(chi_e,pf,eo_conf,appr->poles,appr->nterms,1000000,resiue,residue,0);
  
  //summ all the terms performing appropriate elaboration
  communicate_ev_color_borders(chi_e);
  for(int iterm=0;iterm<appr->nterm;iterm++)
    apply_st2Doe(v_0[iterm],eo+conf,chi_e[iterm]);
  
  //conclude the calculation of the fermionic force
  memset(F,0,loc_vol*sizeof(quad_su3));
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      
    }
  
  //free the back field
  nissa_free(v_o);
  nissa_free(chi_e);
}
