#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "../../inverters/staggered/cg_invert_stD2ee_m2.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/rat_approx.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/mpi.h"

#include "../backfield.h"

#include "rat_expansion_database.cpp"

//compute the norm2 of an even color vector
double eo_color_norm2(color *v)
{
  double norm2=0;
  
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	norm2+=v[ivol][ic][ri]*v[ivol][ic][ri];
  
  return glb_reduce_double(norm2);
}

//put the passed vector to the required norm
double eo_color_normalize(color *out,color *in,double norm)
{
  double fact=sqrt(norm/eo_color_norm2(in));
  
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]=in[ivol][ic][ri]*fact;
  
  set_borders_invalid(out);
  
  return 1/fact;
}

//Return the maximal eigenvalue of the staggered Dirac operator for the passed quark
//assumes that the passed conf already has stag phases inside it
double max_eigenval(quark_content_type &quark_content,quad_su3 **eo_conf,int niters)
{
  communicate_eo_quad_su3_borders(eo_conf);
  
  double eig_max=0;
  
  color *vec_in=nissa_malloc("vec_in",loc_volh+bord_volh,color);
  color *vec_out=nissa_malloc("vec_out",loc_volh,color);
  color *tmp=nissa_malloc("tmp",loc_volh+bord_volh,color);
  
  //generate the random field
  nissa_loc_volh_loop(ivol)
    color_put_to_gauss(vec_in[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),3);
  set_borders_invalid(vec_in);
  
  verbosity_lv3_master_printf("Init norm: %lg\n",eo_color_normalize(vec_in,vec_in,glb_volh*3));
  
  //apply the vector niter times normalizing at each iter
  for(int iter=0;iter<niters;iter++)
    {
      //inv_stD2ee_m2_cg(vec_out,NULL,eo_conf,quark_content.mass,10000,5,1.e-13,vec_in);
      apply_stD2ee(vec_out,eo_conf,tmp,quark_content.mass,vec_in);
      
      //compute the norm
      eig_max=eo_color_normalize(vec_in,vec_out,glb_volh*3);
      
      verbosity_lv2_master_printf("max_eigen search mass %lg, iter=%d, eig=%lg\n",quark_content.mass,iter,eig_max);
    }
  
  nissa_free(tmp);
  nissa_free(vec_out);
  nissa_free(vec_in);
  
  verbosity_lv2_master_printf("max_eigen mass %lg: %lg\n",quark_content.mass,eig_max);
  
  return eig_max;
}

//scale the rational expansion
//assumes that the conf has already stag phases inside
void rootst_eoimpr_scale_expansions(rat_approx_type *rat_exp_pfgen,rat_approx_type *rat_exp_actio,quad_su3 **eo_conf,theory_pars_type *theory_pars)
{
  //loop over each flav
  for(int iflav=0;iflav<theory_pars->nflavs;iflav++)
    {
      //The expansion power for a n-degenerate set of flav is labelled by n-1
      int irexp=theory_pars->quark_content[iflav].deg-1;
      
      //The expansion for a npf pseudofermions is labelled by npf-1
      int ipf=0;//for the moment, only 1 pf
      
      //Set scale factors
      //add the background field
      add_backfield_to_conf(eo_conf,theory_pars->backfield[iflav]);
      double scale=max_eigenval(theory_pars->quark_content[iflav],eo_conf,50)*1.1;
      if(db_rat_exp_min*scale>=pow(theory_pars->quark_content[iflav].mass,2))
	crash("approx not valid, scaled lower range: %lg, min eig: %lg",db_rat_exp_min*scale,pow(theory_pars->quark_content[iflav].mass,2));
      rem_backfield_from_conf(eo_conf,theory_pars->backfield[iflav]);

      double scale_pfgen=pow(scale,db_rat_exp_pfgen_degr[ipf][irexp]);
      double scale_actio=pow(scale,db_rat_exp_actio_degr[ipf][irexp]);

      //scale the rational approximation to generate pseudo-fermions
      rat_exp_pfgen[iflav].exp_power=db_rat_exp_pfgen_degr[ipf][irexp];
      rat_exp_pfgen[iflav].minimum=db_rat_exp_min*scale;
      rat_exp_pfgen[iflav].maximum=db_rat_exp_max*scale;
      rat_exp_pfgen[iflav].cons=db_rat_exp_pfgen_cons[ipf][irexp]*scale_pfgen;

      //scale the rational approximation to compute action (and force)
      rat_exp_actio[iflav].exp_power=db_rat_exp_actio_degr[ipf][irexp];
      rat_exp_actio[iflav].minimum=db_rat_exp_min*scale;
      rat_exp_actio[iflav].maximum=db_rat_exp_max*scale;
      rat_exp_actio[iflav].cons=db_rat_exp_actio_cons[ipf][irexp]*scale_actio;

      for(int iterm=0;iterm<db_rat_exp_nterms;iterm++)
        {
          rat_exp_pfgen[iflav].poles[iterm]=db_rat_exp_pfgen_pole[ipf][irexp][iterm]*scale+pow(theory_pars->quark_content[iflav].mass,2);
          rat_exp_pfgen[iflav].weights[iterm]=db_rat_exp_pfgen_coef[ipf][irexp][iterm]*scale*scale_pfgen;

          rat_exp_actio[iflav].poles[iterm]=db_rat_exp_actio_pole[ipf][irexp][iterm]*scale+pow(theory_pars->quark_content[iflav].mass,2);
          rat_exp_actio[iflav].weights[iterm]=db_rat_exp_actio_coef[ipf][irexp][iterm]*scale*scale_actio;
	}
      
      if(nissa_verbosity>=3)
	{
	  master_printf_rat_approx(&(rat_exp_pfgen[iflav]));
	  master_printf_rat_approx(&(rat_exp_actio[iflav]));
	}
    }
}
