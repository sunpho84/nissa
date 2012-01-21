#pragma once

//compute the norm2 of an even color vector
double eo_color_norm2(color *v)
{
  double norm2=0;
  
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	norm2+=v[ivol][ic][ri]*v[ivol][ic][ri];
  
  return norm2;
}

//put the passed vector to the required norm
double eo_color_normalize(color *v,double norm)
{
  double ori_norm=sqrt(eo_color_norm2(v));
  double inv_fact=norm/ori_norm;
  
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	v[ivol][ic][ri]*=inv_fact;
  
  return ori_norm;
}

//Return the maximal eigenvalue of the staggered Dirac operator for the passed quark
double max_eigenval(quark_content *pars,quad_su3 **eo_conf,int niters)
{
  communicate_eo_gauge_borders(eo_conf[0],eo_conf[1]);
  
  double eig_max;
  
  color *vec=nissa_malloc("vec",(loc_vol+loc_bord)/2,color);
  color *tmp=nissa_malloc("tmp",(loc_vol+loc_bord)/2,color);
  
  //generate the random field
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    color_put_to_gauss(vec[ivol],&(loc_rnd_gen[ivol]),1);  
  
  //apply the vector niter times normalizing at each iter
  for(int iter=0;iter<niters;iter++)
    {
      communicate_ev_color_borders(vec);
      apply_stD2ee(vec,eo_conf,tmp,pars->mass,vec);
      //compute the norm
      eig_max=eo_color_normalize(vec,1);
      
      master_printf("max_eigen search, iter=%d, eig=%lg\n",iter,eig_max);
    }
  
  nissa_free(vec);
  nissa_free(tmp);
  
  return eig_max;
}

//scale the rational expansion
void scale_expansions(rat_approx *rat_exp_pfgen,quad_su3 **eo_conf,quark_content *quark_pars,quad_u1 ***bf,int nquarks)
{
  //loop over each quark
  for(int iquark=0;iquark<nquarks;iquark++)
    {
      //The expansion power for a n-degenerate set of quark is labelled by n-1
      int irexp=quark_pars[iquark].deg-1;
      
      //The expansion for a npf pseudofermions is labelled by npf-1
      int ipf=0;//for the moment, only 1 pf
      
      //Set scale factors
      //add the background field
      /*
	add_backfield_to_conf(eo_conf,u1b);
      double m=max_eigenval(quark_pars,eo_conf);
      rem_backfield_to_conf(eo_conf,u1b);
      */
      double scale=5.6331760852836039;
      double scale_pfgen=pow(scale,db_rat_exp_pfgen_degr[ipf][irexp]);
      master_printf("%lg %lg\n",scale,scale_pfgen);
      
      //scale the rational approximation to generate pseudo-fermions
      rat_exp_pfgen[iquark].exp_power=db_rat_exp_pfgen_degr[ipf][irexp];
      rat_exp_pfgen[iquark].minimum=db_rat_exp_min*scale;
      rat_exp_pfgen[iquark].maximum=db_rat_exp_max*scale;
      rat_exp_pfgen[iquark].cons=db_rat_exp_pfgen_cons[ipf][irexp]*scale_pfgen;

      for(int iterm=0;iterm<db_rat_exp_nterms;iterm++)
        {
          rat_exp_pfgen[iquark].poles[iterm]=db_rat_exp_pfgen_pole[ipf][irexp][iterm]*scale+pow(quark_pars[iquark].mass,2);
          rat_exp_pfgen[iquark].weights[iterm]=db_rat_exp_pfgen_coef[ipf][irexp][iterm]*scale*scale_pfgen;
	}
      
      master_printf_rat_approx(&(rat_exp_pfgen[iquark]));
    }
}
