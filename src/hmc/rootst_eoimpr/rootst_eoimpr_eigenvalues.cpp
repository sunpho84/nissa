#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "communicate/communicate.hpp"
#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "inverters/staggered/cg_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "hmc/backfield.hpp"

#include "rat_expansion_database.hpp"

namespace nissa
{
  //Return the maximal eigenvalue of the staggered Dirac operator for the passed quark
  //assumes that the passed conf already has stag phases inside it
  THREADABLE_FUNCTION_4ARG(max_eigenval, double*,eig_max, quark_content_t*,quark_content, quad_su3**,eo_conf, int,niters)
  {
    niters=10000;
    GET_THREAD_ID();
    communicate_ev_and_od_quad_su3_borders(eo_conf);
    (*eig_max)=0;
    
    color *vec_in1=nissa_malloc("vec_in",loc_volh+bord_volh,color);
    color *vec_in2=nissa_malloc("vec_in",loc_volh+bord_volh,color);
    color *vec_in3=nissa_malloc("vec_in",loc_volh+bord_volh,color);
    color *vec_in4=nissa_malloc("vec_in",loc_volh+bord_volh,color);
    color *vec_out1=nissa_malloc("vec_out",loc_volh,color);
    color *vec_out2=nissa_malloc("vec_out",loc_volh,color);
    color *vec_out3=nissa_malloc("vec_out",loc_volh,color);
    color *vec_out4=nissa_malloc("vec_out",loc_volh,color);
    color *tmp=nissa_malloc("tmp",loc_volh+bord_volh,color);
    
    //generate the random field
    if(IS_MASTER_THREAD)
      NISSA_LOC_VOLH_LOOP(ivol)
      {
	color_put_to_gauss(vec_in1[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),3);
	color_put_to_gauss(vec_in2[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),3);
	color_put_to_gauss(vec_in3[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),3);
	color_put_to_gauss(vec_in4[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),3);
      }
    set_borders_invalid(vec_in1);
    set_borders_invalid(vec_in2);
    set_borders_invalid(vec_in3);
    set_borders_invalid(vec_in4);
    
    //perform initial normalization
    double init_norm1,init_norm2,init_norm3,init_norm4;
    double_vector_normalize(&init_norm1,(double*)vec_in1,(double*)vec_in1,glb_volh*3,6*loc_volh);
    double_vector_normalize(&init_norm2,(double*)vec_in2,(double*)vec_in2,glb_volh*3,6*loc_volh);
    double_vector_normalize(&init_norm3,(double*)vec_in3,(double*)vec_in3,glb_volh*3,6*loc_volh);
    double_vector_normalize(&init_norm4,(double*)vec_in4,(double*)vec_in4,glb_volh*3,6*loc_volh);
    verbosity_lv3_master_printf("Init norm1: %lg\n",init_norm1);
    verbosity_lv3_master_printf("Init norm2: %lg\n",init_norm2);
    verbosity_lv3_master_printf("Init norm3: %lg\n",init_norm3);
    verbosity_lv3_master_printf("Init norm4: %lg\n",init_norm4);
    
    //apply the vector niter times normalizing at each iter
    double eig_max1,eig_max2,eig_max3,eig_max4;
    for(int iter=0;iter<niters;iter++)
      {
	double pr;

	//21
	double_vector_glb_scalar_prod(&pr,(double*)vec_in2,(double*)vec_in1,loc_volh*6);
	double_vector_summ_double_vector_prod_double((double*)vec_in2,(double*)vec_in2,(double*)vec_in1,-pr/(glb_volh*3),loc_volh*6);
	double_vector_normalize(&init_norm1,(double*)vec_in2,(double*)vec_in2,glb_volh*3,6*loc_volh);
	//31
	double_vector_glb_scalar_prod(&pr,(double*)vec_in3,(double*)vec_in1,loc_volh*6);
	double_vector_summ_double_vector_prod_double((double*)vec_in3,(double*)vec_in3,(double*)vec_in1,-pr/(glb_volh*3),loc_volh*6);
	double_vector_normalize(&init_norm1,(double*)vec_in3,(double*)vec_in3,glb_volh*3,6*loc_volh);
	//32
	double_vector_glb_scalar_prod(&pr,(double*)vec_in3,(double*)vec_in2,loc_volh*6);
	double_vector_summ_double_vector_prod_double((double*)vec_in3,(double*)vec_in3,(double*)vec_in2,-pr/(glb_volh*3),loc_volh*6);
	double_vector_normalize(&init_norm1,(double*)vec_in3,(double*)vec_in3,glb_volh*3,6*loc_volh);
	//41
	double_vector_glb_scalar_prod(&pr,(double*)vec_in4,(double*)vec_in1,loc_volh*6);
	double_vector_summ_double_vector_prod_double((double*)vec_in4,(double*)vec_in4,(double*)vec_in1,-pr/(glb_volh*3),loc_volh*6);
	double_vector_normalize(&init_norm1,(double*)vec_in4,(double*)vec_in4,glb_volh*3,6*loc_volh);
	//42
	double_vector_glb_scalar_prod(&pr,(double*)vec_in4,(double*)vec_in2,loc_volh*6);
	double_vector_summ_double_vector_prod_double((double*)vec_in4,(double*)vec_in4,(double*)vec_in2,-pr/(glb_volh*3),loc_volh*6);
	double_vector_normalize(&init_norm1,(double*)vec_in4,(double*)vec_in4,glb_volh*3,6*loc_volh);
	//43
	double_vector_glb_scalar_prod(&pr,(double*)vec_in4,(double*)vec_in3,loc_volh*6);
	double_vector_summ_double_vector_prod_double((double*)vec_in4,(double*)vec_in4,(double*)vec_in3,-pr/(glb_volh*3),loc_volh*6);
	double_vector_normalize(&init_norm1,(double*)vec_in4,(double*)vec_in4,glb_volh*3,6*loc_volh);

	//inv_stD2ee_m2_cg(vec_out,NULL,eo_conf,sqr(quark_content.mass),10000,5,1.e-13,vec_in);
	apply_stD2ee_m2(vec_out1,eo_conf,tmp,quark_content->mass*quark_content->mass,vec_in1);
	apply_stD2ee_m2(vec_out2,eo_conf,tmp,quark_content->mass*quark_content->mass,vec_in2);
	apply_stD2ee_m2(vec_out3,eo_conf,tmp,quark_content->mass*quark_content->mass,vec_in3);
	apply_stD2ee_m2(vec_out4,eo_conf,tmp,quark_content->mass*quark_content->mass,vec_in4);
	
	//compute the norm
	double_vector_normalize(&eig_max1,(double*)vec_in1,(double*)vec_out1,glb_volh*3,6*loc_volh);
	double_vector_normalize(&eig_max2,(double*)vec_in2,(double*)vec_out2,glb_volh*3,6*loc_volh);
	double_vector_normalize(&eig_max3,(double*)vec_in3,(double*)vec_out3,glb_volh*3,6*loc_volh);
	double_vector_normalize(&eig_max4,(double*)vec_in4,(double*)vec_out4,glb_volh*3,6*loc_volh);
	//verbosity_lv2_
	master_printf("max_eigen search mass %lg, iter=%d, eig1=%16.16lg, eig2=%16.16lg, eig3=%16.16lg, eig4=%16.16lg\n",quark_content->mass,iter,eig_max1,eig_max2,eig_max3,eig_max4);
      }
    
    nissa_free(tmp);
    nissa_free(vec_out1);
    nissa_free(vec_out2);
    nissa_free(vec_out3);
    nissa_free(vec_out4);
    nissa_free(vec_in1);
    nissa_free(vec_in2);
    nissa_free(vec_in3);
    nissa_free(vec_in4);
    
    (*eig_max)=std::max(eig_max1,eig_max2);
    
    verbosity_lv2_master_printf("max_eigen mass %lg: %16.16lg %16.16g\n",quark_content->mass,eig_max1,eig_max2);
  }
  THREADABLE_FUNCTION_END

  //scale the rational expansion
  //assumes that the conf has already stag phases inside
  THREADABLE_FUNCTION_5ARG(rootst_eoimpr_scale_expansions, rat_approx_t*,rat_exp_pfgen, rat_approx_t*,rat_exp_actio, quad_su3**,eo_conf, theory_pars_t*,theory_pars, int*,npfs)
  {
    GET_THREAD_ID();
    
    //loop over each flav
    for(int iflav=0;iflav<theory_pars->nflavs;iflav++)
      {
	//The expansion power for a n-degenerate set of flav is labelled by n-1
	int irexp=theory_pars->quark_content[iflav].deg-1;
	
	//The expansion for a npf pseudofermions is labelled by npf-1
	int ipf=npfs[iflav]-1;
	
	//Set scale factors
	//add the background field
	add_backfield_to_conf(eo_conf,theory_pars->backfield[iflav]);
	double scale;
	max_eigenval(&scale,&(theory_pars->quark_content[iflav]),eo_conf,50);
	scale*=1.1;
	if(db_rat_exp_min*scale>=pow(theory_pars->quark_content[iflav].mass,2))
	  crash("approx not valid, scaled lower range: %lg, min eig: %lg",db_rat_exp_min*scale,pow(theory_pars->quark_content[iflav].mass,2));
	rem_backfield_from_conf(eo_conf,theory_pars->backfield[iflav]);
	
	if(IS_MASTER_THREAD)
	  {
	    double scale_pfgen=pow(scale,db_rat_exp_pfgen_degr[ipf][irexp]);
	    double scale_actio=pow(scale,db_rat_exp_actio_degr[ipf][irexp]);
	    
	    //scale the rational approximation to generate pseudo-fermions
	    rat_exp_pfgen[iflav].minimum=db_rat_exp_min*scale;
	    rat_exp_pfgen[iflav].maximum=db_rat_exp_max*scale;
	    rat_exp_pfgen[iflav].cons=db_rat_exp_pfgen_cons[ipf][irexp]*scale_pfgen;
	    
	    //scale the rational approximation to compute action (and force)
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
	  }
	
	//sync
	THREAD_BARRIER();
	
	if(VERBOSITY_LV3)
	  {
	    master_printf_rat_approx(&(rat_exp_pfgen[iflav]));
	    master_printf_rat_approx(&(rat_exp_actio[iflav]));
	  }
      }
  }
  THREADABLE_FUNCTION_END
  
}
