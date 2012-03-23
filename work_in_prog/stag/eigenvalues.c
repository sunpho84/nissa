#pragma once

//compute the norm2 of an even color vector
double eo_color_norm2(color *v)
{
  double loc_norm2=0;
  
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	loc_norm2+=v[ivol][ic][ri]*v[ivol][ic][ri];
  
  double norm2;
  MPI_Allreduce(&loc_norm2,&norm2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return norm2;
}

//put the passed vector to the required norm
double eo_color_normalize(color *out,color *in,double norm)
{
  double fact=sqrt(norm/eo_color_norm2(in));
  
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]=in[ivol][ic][ri]*fact;
  
  return 1/fact;
}

//Return the maximal eigenvalue of the staggered Dirac operator for the passed quark
//assumes that the passed conf already has stag phases inside it
double max_eigenval(quark_content *pars,quad_su3 **eo_conf,int niters)
{
  communicate_eo_quad_su3_borders(eo_conf[0],eo_conf[1]);
  
  double eig_max;
  
  color *vec_in=nissa_malloc("vec_in",loc_volh+loc_bordh,color);
  color *vec_out=nissa_malloc("vec_out",loc_volh+loc_bordh,color);
  color *tmp=nissa_malloc("tmp",loc_volh+loc_bordh,color);
  
  //generate the random field
  nissa_loc_volh_loop(ivol)
    color_put_to_gauss(vec_in[ivol],&(loc_rnd_gen[loclx_of_loceo[0][ivol]]),3);
  
  /*
  //debug
  static int a=0;
  master_printf("Debug, loading initial eig\n");
  if(a==0) read_e_color(vec_in,"dat/Eig_init1");
  else     read_e_color(vec_in,"dat/Eig_init2");
  a++;
  */
  
  //apply the vector niter times normalizing at each iter
  for(int iter=0;iter<niters;iter++)
    {
      int n=loceo_of_loclx[loclx_of_coord_list(0,7,0,0)];
      communicate_ev_color_borders(vec_in);
      n=loceo_neighup[1][n][1];
      
      int r,nx,neo;
      coords g={0,8,0,0};
      get_loclx_and_rank_of_coord(&nx,&r,g);
      neo=loceo_of_loclx[nx];
      
      apply_stD2ee(vec_out,eo_conf,tmp,pars->mass,vec_in);
      
      //compute the norm
      eig_max=eo_color_normalize(vec_in,vec_out,glb_volh*3);
      
      master_printf("max_eigen search, iter=%d, eig=%lg\n",iter,eig_max);
    }
  
  nissa_free(tmp);
  nissa_free(vec_out);
  nissa_free(vec_in);
  
  return eig_max;
}

//scale the rational expansion
//assumes that the conf has already stag phases inside
void scale_expansions(rat_approx *rat_exp_pfgen,rat_approx *rat_exp_actio,quad_su3 **eo_conf,quark_content *quark_pars,quad_u1 ***bf,int nquarks)
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
      add_backfield_to_conf(eo_conf,bf[iquark]);
      double scale=max_eigenval(quark_pars,eo_conf,50)*1.1;
      rem_backfield_from_conf(eo_conf,bf[iquark]);

      /*
      ///// DEBUG //////
      double exp_scale=5.6331760852836039;
      master_printf("Debug: scaling with %16.16lg, expected: %16.16lg\n",scale,exp_scale);
      master_printf("Forcing to %16.16lg\n",exp_scale);
      scale=exp_scale;
      */
      
      double scale_pfgen=pow(scale,db_rat_exp_pfgen_degr[ipf][irexp]);
      double scale_actio=pow(scale,db_rat_exp_actio_degr[ipf][irexp]);

      //scale the rational approximation to generate pseudo-fermions
      rat_exp_pfgen[iquark].exp_power=db_rat_exp_pfgen_degr[ipf][irexp];
      rat_exp_pfgen[iquark].minimum=db_rat_exp_min*scale;
      rat_exp_pfgen[iquark].maximum=db_rat_exp_max*scale;
      rat_exp_pfgen[iquark].cons=db_rat_exp_pfgen_cons[ipf][irexp]*scale_pfgen;

      //scale the rational approximation to compute action (and force)
      rat_exp_actio[iquark].exp_power=db_rat_exp_actio_degr[ipf][irexp];
      rat_exp_actio[iquark].minimum=db_rat_exp_min*scale;
      rat_exp_actio[iquark].maximum=db_rat_exp_max*scale;
      rat_exp_actio[iquark].cons=db_rat_exp_actio_cons[ipf][irexp]*scale_actio;

      for(int iterm=0;iterm<db_rat_exp_nterms;iterm++)
        {
          rat_exp_pfgen[iquark].poles[iterm]=db_rat_exp_pfgen_pole[ipf][irexp][iterm]*scale+pow(quark_pars[iquark].mass,2);
          rat_exp_pfgen[iquark].weights[iterm]=db_rat_exp_pfgen_coef[ipf][irexp][iterm]*scale*scale_pfgen;

          rat_exp_actio[iquark].poles[iterm]=db_rat_exp_actio_pole[ipf][irexp][iterm]*scale+pow(quark_pars[iquark].mass,2);
          rat_exp_actio[iquark].weights[iterm]=db_rat_exp_actio_coef[ipf][irexp][iterm]*scale*scale_actio;
	}
      
      master_printf_rat_approx(&(rat_exp_pfgen[iquark]));
      master_printf_rat_approx(&(rat_exp_actio[iquark]));
    }
}
