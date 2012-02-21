#pragma once

//Compute the fermionic force the rooted staggered e/o improved theory.
//Passed conf must NOT contain the backfield.
//Of the result still need to be taken the TA
//The approximation need to be already scaled, and must contain physical mass term
void summ_the_rootst_eoimpr_force(quad_su3 **F,quad_su3 **eo_conf,quad_u1 **u1b,color *pf,rat_approx *appr,double residue)
{
  master_printf("Computing quark force\n");
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

  //debug
  static int b=0;
  
  //invert the various terms
  master_printf_rat_approx(appr);
  inv_stD2ee_cgmm2s(chi_e,pf,eo_conf,appr->poles,appr->nterms,1000000,residue,residue,0);
  
  //summ all the terms performing appropriate elaboration
  //possible improvement by communicating more borders together
  for(int iterm=0;iterm<appr->nterms;iterm++)
    {
      communicate_ev_color_borders(chi_e[iterm]);
      apply_stDoe(v_o[iterm],eo_conf,chi_e[iterm]);
    }
  
  //debug
  master_printf("Debug,reading the first term\n");
  color *temp=nissa_malloc("temp",loc_volh,color);
  if(b==0) read_e_color(temp,"dat/quark_force_term_1_1");
  else     read_e_color(temp,"dat/quark_force_term_2_1");
  double n2=0;
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  double a=temp[ivol][ic][ri]-chi_e[0][ivol][ic][ri];
	  n2+=a*a;
	  master_printf("%d %d %d  %lg %lg\n",ivol,ic,ri,temp[ivol][ic][ri],chi_e[0][ivol][ic][ri]);
	}
  n2=sqrt(n2/(loc_volh*3));
  master_printf("Diff of first term: %lg\n",n2);
  if(n2>=1.e-4) crash("norm error");
  
  master_printf("Debug,reading the first vterm\n");
  if(b==0) read_o_color(temp,"dat/quark_force_vterm_1_1");
  else     read_o_color(temp,"dat/quark_force_vterm_2_1");
  b++;
  n2=0;
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  double a=temp[ivol][ic][ri]-v_o[0][ivol][ic][ri];
	  n2+=a*a;
	  master_printf("%d %d %d  %lg %lg\n",ivol,ic,ri,temp[ivol][ic][ri],v_o[0][ivol][ic][ri]);
	}
  n2=sqrt(n2/(loc_volh*3));
  master_printf("Diff of first term: %lg\n",n2);
  if(n2>=1.e-4) crash("norm error");
  
  //remove the background fields
  rem_backfield_from_conf(eo_conf,u1b);
  
  //conclude the calculation of the fermionic force
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int mu=0;mu<4;mu++)
      for(int iterm=0;iterm<appr->nterms;iterm++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    {
	      //chpot to be added
	      complex temp;
	      
	      //this is for ivol=EVN
	      unsafe_complex_conj2_prod(temp,v_o[iterm][loceo_neighup[EVN][ivol][mu]][ic1],chi_e[iterm][ivol][ic2]);
	      complex_summ_the_prod_double(F[EVN][ivol][mu][ic1][ic2],temp,appr->weights[iterm]);
	      
	      //this is for ivol=ODD
	      unsafe_complex_conj2_prod(temp,chi_e[iterm][loceo_neighup[ODD][ivol][mu]][ic1],v_o[iterm][ivol][ic2]);
	      complex_subt_the_prod_double(F[ODD][ivol][mu][ic1][ic2],temp,appr->weights[iterm]);
	    }
  
  //free the back field
  nissa_free(v_o[0]);
  nissa_free(chi_e[0]);
}

//Compute the full force for the rooted staggered theory with nfl flavors
//Conf and pf borders are assumed to have been already communicated
void full_rootst_eoimpr_force(quad_su3 **F,quad_su3 **conf,double beta,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue)
{
  //First of all compute gluonic part of the force
  add_backfield_to_conf(conf,u1b[0]);
  wilson_force(F,conf,beta);
  rem_backfield_from_conf(conf,u1b[0]);
  
  //Then summ the contributes coming from each flavor
  for(int ifl=0;ifl<nfl;ifl++)
    summ_the_rootst_eoimpr_force(F,conf,u1b[ifl],pf[ifl],&(appr[ifl]),residue);
  
  //debug
  quad_su3 *temp[2];
  temp[0]=nissa_malloc("temp",loc_vol,quad_su3);
  temp[1]=temp[0]+loc_volh;
  read_ildg_gauge_conf_and_split_into_eo_parts(temp,"dat/total_force_pre_ta");
  double n2=0;
  for(int eo=0;eo<2;eo++)
    for(int ivol=0;ivol<loc_volh;ivol++)
      for(int mu=0;mu<4;mu++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    for(int ri=0;ri<2;ri++)
	      {
		double a=temp[eo][ivol][mu][ic1][ic2][ri]-F[eo][ivol][mu][ic1][ic2][ri];
		master_printf("%d %d %d %d %d %d  %lg %lg\n",eo,ivol,mu,ic1,ic2,ri,temp[eo][ivol][mu][ic1][ic2][ri],F[eo][ivol][mu][ic1][ic2][ri]);
		n2+=a*a;
	      }
  n2/=loc_vol*4*9;
  n2=sqrt(n2);
  master_printf("Total force norm diff: %lg\n",n2);
  if(n2>1.e-16) crash("quark and wilson force computatin failed");
  
  //Finish the computation multiplying for the conf and aking TA
  for(int eo=0;eo<2;eo++)
    for(int ivol=0;ivol<loc_volh;ivol++)
      for(int mu=0;mu<4;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[eo][ivol][mu],F[eo][ivol][mu]);
	  su3_traceless_anti_hermitian_part(F[eo][ivol][mu],temp);
	}
}
