#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../base/random.h"
#include "../../base/thread_macros.h"
#include "../../base/vectors.h"
#include "../../communicate/communicate.h"
#include "../../geometry/geometry_eo.h"
#include "../../hmc/backfield.h"
#include "../../new_types/complex.h"
#include "../../inverters/staggered/cg_invert_stD.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/mpi_routines.h"
#ifdef USE_THREADS
 #include "../../routines/thread.h"
#endif

//get a propagator
THREADABLE_FUNCTION_6ARG(get_propagator, color**,prop, quad_su3**,conf, quad_u1**,u1b, double,m, double,residue, color**,source)
{
  addrem_stagphases_to_eo_conf(conf);
  add_backfield_to_conf(conf,u1b);
  inv_stD_cg(prop,conf,m,10000,5,residue,source);
  rem_backfield_from_conf(conf,u1b);
  addrem_stagphases_to_eo_conf(conf);
}}

//compute the chiral condensate
THREADABLE_FUNCTION_5ARG(chiral_condensate, double*,cond, quad_su3**,conf, quad_u1**,u1b, quark_content_t*,quark, double,residue)
{
  //allocate
  color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
  color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
  
  //generate the source and the propagator
  generate_fully_undiluted_eo_source(rnd,RND_GAUSS,-1);
  get_propagator(chi,conf,u1b,quark->mass,residue,rnd);
  
  //summ the scalar prod of EVN and ODD parts
  double temp[2];
  for(int eo=0;eo<2;eo++)
    double_vector_glb_scalar_prod(temp+eo,(double*)(rnd[eo]),(double*)(chi[eo]),3*loc_volh*2);
  
  //add normalization: deg/4vol
  (*cond)=(temp[EVN]+temp[ODD])*quark->deg/(4.0*glb_vol);
  
  //free
  for(int par=0;par<2;par++)
    {
      nissa_free(rnd[par]);
      nissa_free(chi[par]);
    }
}}

//measure chiral cond
void measure_chiral_cond(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created)
{
  FILE *file=open_file(theory_pars.chiral_cond_pars.path,conf_created?"w":"a");

  master_fprintf(file,"%d",iconf);
  
  //measure the condensate for each quark
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    {
      double cond=0;
      
      //loop over hits
      int nhits=theory_pars.chiral_cond_pars.nhits;
      for(int hit=0;hit<nhits;hit++)
        {
          verbosity_lv1_master_printf("Evaluating chiral condensate for flavor %d/%d, nhits %d/%d\n",iflav+1,theory_pars.nflavs,hit+1,nhits);
          
          //compute and summ
          double temp;
          chiral_condensate(&temp,conf,theory_pars.backfield[iflav],theory_pars.quark_content+iflav,theory_pars.chiral_cond_pars.residue);
	  cond+=temp/nhits;
        }
      
      master_fprintf(file,"\t%+016.16lg",cond);
    }

  master_fprintf(file,"\n");
  
  if(rank==0) fclose(file);
}

//compute the magnetization
THREADABLE_FUNCTION_5ARG(magnetization, complex*,magn, quad_su3**,conf, quad_u1**,u1b, quark_content_t*,quark, double,residue)
{
  GET_THREAD_ID();
  
  //fixed to Z magnetization
  int mu=1,nu=2;
  
  //allocate source and propagator
  color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
  color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
  
  //we need to store phases
  coords *arg=nissa_malloc("arg",loc_vol+bord_vol,coords);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol+bord_vol)
    get_args_of_one_over_L2_quantization(arg[ivol],ivol,mu,nu);
  
  //array to store magnetization on single site (actually storing backward contrib at displaced site)
  complex *point_magn=nissa_malloc("app",loc_vol,complex);
  vector_reset(point_magn);
  
  //generate the source and the propagator
  generate_fully_undiluted_eo_source(rnd,RND_GAUSS,-1);
  
  //we add stagphases and backfield externally because we need them for derivative
  addrem_stagphases_to_eo_conf(conf);
  add_backfield_to_conf(conf,u1b);
  
  //invert
  inv_stD_cg(chi,conf,quark->mass,10000,5,residue,rnd);
  communicate_ev_and_od_color_borders(chi);
  
  //summ the results of the derivative
  for(int par=0;par<2;par++)
    NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
      {
	int ivol=loclx_of_loceo[par][ieo];
	
	//summ the contribution of the derivative in mu and nu directions
	int rho_list[2]={mu,nu};
	for(int irho=0;irho<2;irho++)
	  {
	    int rho=rho_list[irho];
	    
	    int iup_eo=loceo_neighup[par][ieo][rho];
	    int idw_eo=loceo_neighdw[par][ieo][rho];
	    int idw_lx=loclx_neighdw[ivol][rho];
	    
	    color v;
	    complex t;
	    
	    //forward derivative
	    unsafe_su3_prod_color(v,conf[par][ieo][rho],chi[!par][iup_eo]);
	    color_scalar_prod(t,v,rnd[par][ieo]);
	    complex_summ_the_prod_double(point_magn[ivol],t,arg[ivol][rho]);
	    
	    //backward derivative: note that we should multiply for -arg*(-U^+)
	    unsafe_su3_dag_prod_color(v,conf[!par][idw_eo][rho],chi[!par][idw_eo]);
	    color_scalar_prod(t,v,rnd[par][ieo]);
	    complex_summ_the_prod_double(point_magn[ivol],t,arg[idw_lx][rho]);
	  }
      }
  
  //remove stag phases and u1 field, and automatically barrier before collapsing
  rem_backfield_from_conf(conf,u1b);
  addrem_stagphases_to_eo_conf(conf);
  
  //reduce across all nodes and threads
  complex temp;
  complex_vector_glb_collapse(temp,point_magn,loc_vol);
  
  //add normalization, corresponding to all factors relative to derivative with respects to "b": 
  //-quark_deg/4 coming from the determinant
  //-1/vol coming from stochastic trace
  //-1/2 coming from dirac operator
  //-i*2*quark_charge*M_PI/glb_size[mu]/glb_size[nu] coming EM potential prefactor in front of "b"
  //and a minus because F=-logZ
  unsafe_complex_prod_idouble(*magn,temp,-quark->deg*2*M_PI*quark->charge/(4.0*glb_vol*2*glb_size[mu]*glb_size[nu]));
  
  //free
  for(int par=0;par<2;par++)
    {
      nissa_free(rnd[par]);
      nissa_free(chi[par]);
    }
  nissa_free(point_magn);
  nissa_free(arg);
}}

//measure magnetization
void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created)
{
  FILE *file=open_file(theory_pars.magnetization_pars.path,conf_created?"w":"a");

  master_fprintf(file,"%d",iconf);
  
  //measure magnetization for each quark
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    {
      complex magn={0,0};
      
      //loop over hits
      int nhits=theory_pars.magnetization_pars.nhits;
      for(int hit=0;hit<nhits;hit++)
        {
          verbosity_lv1_master_printf("Evaluating magnetization for flavor %d/%d, nhits %d/%d\n",iflav+1,theory_pars.nflavs,hit+1,nhits);
          
          //compute and summ
          complex temp;
          magnetization(&temp,conf,theory_pars.backfield[iflav],theory_pars.quark_content+iflav,
			theory_pars.magnetization_pars.residue);
          complex_summ_the_prod_double(magn,temp,1.0/nhits);
        }
      
      master_fprintf(file,"\t%+016.16lg \t%+016.16lg",magn[RE],magn[IM]);
    }

  master_fprintf(file,"\n");
  
  if(rank==0) fclose(file);
}
 
//compute the local pseudoscalar correlator
void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created)
{
  FILE *file=open_file(theory_pars.pseudo_corr_pars.path,conf_created?"w":"a");
  
  int nflavs=theory_pars.nflavs;
  
  //allocate source and propagators
  color *source[2]={nissa_malloc("source_e",loc_volh+bord_volh,color),nissa_malloc("source_o",loc_volh+bord_volh,color)};
  color *prop[nflavs][2];
  for(int iflav=0;iflav<nflavs;iflav++)
    for(int EO=0;EO<2;EO++)
      prop[iflav][EO]=nissa_malloc("prop",loc_volh+bord_volh,color);
  
  //allocate local and global contraction
  complex *loc_contr=nissa_malloc("loc_contr",glb_size[0]*nflavs*(nflavs+1)/2,complex);
  complex *glb_contr=nissa_malloc("glb_contr",glb_size[0]*nflavs*(nflavs+1)/2,complex);
  vector_reset(loc_contr);
  
  //loop over the hits
  int nhits=theory_pars.pseudo_corr_pars.nhits;
  for(int hit=0;hit<nhits;hit++)
    {
      //generate the source on an even site
      int twall=(int)rnd_get_unif(&glb_rnd_gen,0,glb_size[0]/2)*2;
      generate_fully_undiluted_eo_source(source,RND_Z4,twall);
      //filter_hypercube_origin_sites(source); //this is in conjunction with factor "8"
      
      //compute propagators
      for(int iflav=0;iflav<nflavs;iflav++)
	get_propagator(prop[iflav],conf,theory_pars.backfield[iflav],theory_pars.quark_content[iflav].mass,
		       theory_pars.pseudo_corr_pars.residue,source);
      
      //contract the propagators
      int icombo=0;
      for(int iflav=0;iflav<nflavs;iflav++)
	for(int jflav=0;jflav<=iflav;jflav++)
	  {
	    for(int eo=0;eo<2;eo++)
	      nissa_loc_volh_loop(ieo)
	      {
		int ilx=loclx_of_loceo[eo][ieo];
		int t=(glb_coord_of_loclx[ilx][0]+glb_size[0]-twall)%glb_size[0];
		for(int ic=0;ic<3;ic++)
		  complex_summ_the_conj2_prod(loc_contr[icombo*glb_size[0]+t],prop[iflav][eo][ieo][ic],prop[jflav][eo][ieo][ic]);
	      }
	    icombo++;
	  }
    }

  //reduce
  glb_reduce_complex_vect(glb_contr,loc_contr,glb_size[0]*nflavs*(nflavs+1)/2);
  
  //print
  double norm=nhits*glb_vol/(/*8**/glb_size[0]);
  int icombo=0;
  for(int iflav=0;iflav<nflavs;iflav++)
    for(int jflav=0;jflav<=iflav;jflav++)
      {
	master_fprintf(file," # iconf %d , m1 = %lg , m2 = %lg\n",iconf,theory_pars.quark_content[iflav].mass,theory_pars.quark_content[jflav].mass);
	
        for(int t=0;t<glb_size[0];t++)
          master_fprintf(file,"%d %+016.16lg\n",t,glb_contr[icombo*glb_size[0]+t][RE]/norm);
        icombo++;
      }
  
  //free everything
  for(int EO=0;EO<2;EO++)
    {    
      for(int iflav=0;iflav<nflavs;iflav++)
        nissa_free(prop[iflav][EO]);
      nissa_free(source[EO]);
    }
  nissa_free(loc_contr);
  nissa_free(glb_contr);
  
  if(rank==0) fclose(file);
}

/*
//compute the local pseudoscalar correlator
THREADABLE_FUNCTION_4ARG(measure_time_meson_corr, quad_su3**,conf, theory_pars_t*,theory_pars, int,iconf, int,conf_created)
{
  GET_THREAD_ID();
  
  int nflavs=theory_pars->nflavs;
  
  //allocate source and propagators
  color *ori_source[2]={nissa_malloc("ori_source_e",loc_volh,color),nissa_malloc("ori_source_o",loc_volh,color)};
  color *source[2]={nissa_malloc("source_e",loc_volh+bord_volh,color),nissa_malloc("source_o",loc_volh+bord_volh,color)};
  color *prop[16][nflavs][2];
  for(int icube=0;icube<16;icube++)
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int EO=0;EO<2;EO++)
	prop[icube][iflav][EO]=nissa_malloc("prop",loc_volh+bord_volh,color);
  
  //generate the source on hypercube vertex
  int twall=(int)rnd_get_unif(&glb_rnd_gen,0,glb_size[0]/2)*2;
  generate_fully_undiluted_eo_source(ori_source,RND_Z4,twall);
  filter_hypercube_origin_sites(ori_source);
      
  //compute propagators
  for(int icube=0;icube<16;icube++)
    {
      //transport the source
      vector_reset(source[EVN]);vector_reset(source[ODD]);      
      NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
	if(hyp_parity(ivol_lx)==icube)
	  {
	    int ipar=loclx_parity[ivol_lx];
	    int ivol_eo=loceo_of_loclx[ivol_lx];
	    
	    int hyp_vertex=hyp_vertex_of_loclx(ivol_lx);
	  }
      
      for(int iflav=0;iflav<nflavs;iflav++)
	get_propagator(prop[icube][iflav],conf,theory_pars->backfield[iflav],theory_pars->quark_content[iflav].mass,
		       theory_pars->pseudo_corr_pars.residue,source);
    }
  
  //free everything  
  for(int icube=0;icube<16;icube++)
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int EO=0;EO<2;EO++)
	nissa_free(prop[icube][iflav][EO]);
  nissa_free(source[EVN]);nissa_free(source[ODD]);
  nissa_free(ori_source[EVN]);nissa_free(ori_source[ODD]);
  nissa_free(shortest[EVN]);nissa_free(shortest[ODD]); 
}}
*/
