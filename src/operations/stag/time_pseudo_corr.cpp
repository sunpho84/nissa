#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //compute the local pseudoscalar correlator in "time" direction (that can be all but time)
  void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,pseudo_corr_meas_pars_t &meas_pars,int iconf,int conf_created,int dir)
  {
    char dir_name[5]="txyz";
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    int nflavs=theory_pars.nflavs;
    
    //allocate source
    color *source[2]={nissa_malloc("source_e",loc_volh+bord_volh,color),
		      nissa_malloc("source_o",loc_volh+bord_volh,color)};
    
    //allocate sink
    color *prop[nflavs][2];
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int EO=0;EO<2;EO++)
	prop[iflav][EO]=nissa_malloc("prop",loc_volh+bord_volh,color);
    
    //allocate local and global contraction
    complex *loc_contr=nissa_malloc("loc_contr",glb_size[dir]*nflavs*(nflavs+1)/2,complex);
    complex *glb_contr=nissa_malloc("glb_contr",glb_size[dir]*nflavs*(nflavs+1)/2,complex);
    vector_reset(loc_contr);
    
    //loop over the hits
    int nhits=meas_pars.nhits;
    for(int hit=0;hit<nhits;hit++)
      {
	verbosity_lv2_master_printf("Evaluating pseudoscalar %c correlator, hit %d/%d\n",dir_name[dir],hit+1,nhits);
	
	//generate the source on an even site
	int twall=(int)rnd_get_unif(&glb_rnd_gen,0,glb_size[dir]/2)*2;
	generate_fully_undiluted_eo_source(source,RND_Z4,twall,dir);
	//filter_hypercube_origin_sites(source); //this is in conjunction with factor "8"
	
	//compute M^-1
	for(int iflav=0;iflav<nflavs;iflav++)
	  mult_Minv(prop[iflav],conf,&theory_pars,iflav,meas_pars.residue,source);
	
	//contract
	int icombo=0;
	for(int iflav=0;iflav<nflavs;iflav++)
	  for(int jflav=0;jflav<=iflav;jflav++)
	    {
	      for(int eo=0;eo<2;eo++)
		NISSA_LOC_VOLH_LOOP(ieo)
		{
		  int ilx=loclx_of_loceo[eo][ieo];
		  int t=(glb_coord_of_loclx[ilx][dir]+glb_size[dir]-twall)%glb_size[dir];
		  for(int ic=0;ic<3;ic++)
		    complex_summ_the_conj2_prod(loc_contr[icombo*glb_size[dir]+t],
						prop[iflav][eo][ieo][ic],prop[jflav][eo][ieo][ic]);
		}
	      icombo++;
	    }
      }
    
    //reduce
    glb_nodes_reduce_complex_vect(glb_contr,loc_contr,glb_size[dir]*nflavs*(nflavs+1)/2);
    
    //print
    double norm=nhits*glb_vol/(/*8**/glb_size[dir]);
    int icombo=0;
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int jflav=0;jflav<=iflav;jflav++)
	{
	  master_fprintf(file," # conf %d , \'%c\' dir ; flv1 = %d , m1 = %lg ; flv2 = %d , m2 = %lg\n",
			 iconf,dir_name[dir],
			 iflav,theory_pars.quark_content[iflav].mass,jflav,theory_pars.quark_content[jflav].mass);
	  
	  for(int t=0;t<glb_size[dir];t++)
	    master_fprintf(file,"%d %+016.16lg\n",t,glb_contr[icombo*glb_size[dir]+t][RE]/norm);
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
  
  //allocate source and sink
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
      
  //compute M^-1
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
      mult_Minv(prop[icube][iflav],conf,theory_pars->backfield[iflav],theory_pars->quark_content[iflav].mass,
      theory_pars->pseudo_corr_meas_pars.residue,source);
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
}
