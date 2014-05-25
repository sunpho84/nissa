#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/math_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //get a propagator
  THREADABLE_FUNCTION_6ARG(get_propagator, color**,prop, quad_su3**,conf, quad_u1**,u1b, double,m, double,residue, color**,source)
  {
    addrem_stagphases_to_eo_conf(conf);
    add_backfield_to_conf(conf,u1b);
    do inv_stD_cg(prop,conf,m,100000,residue,source);
    while(issued_cg_warning);
    rem_backfield_from_conf(conf,u1b);
    addrem_stagphases_to_eo_conf(conf);
  }
  THREADABLE_FUNCTION_END
  
  //holds the putpourri in a clean way
  struct fermionic_putpourri_t
  {
    complex chiral_cond;
    complex chiral_cond_susc;
    complex energy_dens;
    complex quark_dens;
    complex quark_dens_susc;
    complex pressure_dens;
    void reset()
    {
      for(int ri=0;ri<2;ri++)
	chiral_cond[ri]=
	  chiral_cond_susc[ri]=
	  energy_dens[ri]=
	  quark_dens[ri]=
	  quark_dens_susc[ri]=
	  pressure_dens[ri]=0;
    }
    fermionic_putpourri_t() {reset();}
  };

  //compute the matrix element of the derivative of the dirac operator between two vectors
  //forward and backward derivative are stored separately, for a reason
  void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result)
  {
    GET_THREAD_ID();
    
    color **right_fw_bw[2]={right,left};

    for(int fw_bw=0;fw_bw<2;fw_bw++)
      {
	vector_reset(point_result);
	for(int par=0;par<2;par++)
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      color v;
	      unsafe_su3_prod_color(v,conf[par][ieo][mu],right_fw_bw[fw_bw][!par][loceo_neighup[par][ieo][mu]]);
	      complex t;
	      if(fw_bw==0) color_scalar_prod(t,v,right_fw_bw[!fw_bw][par][ieo]);
	      else         color_scalar_prod(t,right_fw_bw[!fw_bw][par][ieo],v);
	      complex_summassign(point_result[loclx_of_loceo[par][ieo]],t);
	    }
	THREAD_BARRIER();
	complex_vector_glb_collapse(res_fw_bw[fw_bw],point_result,loc_vol);
      }
  }
    
  //compute the fermionic putpourri for a single conf and hit
  THREADABLE_FUNCTION_6ARG(fermionic_putpourri, fermionic_putpourri_t*,putpourri, quad_su3**,conf, quad_u1**,u1b, quark_content_t*,quark, double,residue, int,comp_susc)
  {
    GET_THREAD_ID();
    
    THREAD_BARRIER();
    
    //allocate
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    color *chi1[2]={nissa_malloc("chi1_EVN",loc_volh+bord_volh,color),nissa_malloc("chi1_ODD",loc_volh+bord_volh,color)};
    color *chi2[2],*app[2];
    if(comp_susc)
      for(int par=0;par<2;par++)
	{
	  chi2[par]=nissa_malloc("chi2_EO",loc_volh+bord_volh,color);
	  app[par]=nissa_malloc("app_EO",loc_volh+bord_volh,color);
	}
    
    //generate the source
    generate_fully_undiluted_eo_source(rnd,RND_GAUSS,-1);
    
    //we add stagphases and backfield externally because we need them for derivative
    addrem_stagphases_to_eo_conf(conf);
    add_backfield_to_conf(conf,u1b);
    
    //invert
    inv_stD_cg(chi1,conf,quark->mass,100000,residue,rnd);
    communicate_ev_and_od_color_borders(chi1);
    
    //invert for chi2 if susceptivity needed
    if(comp_susc)
      {
	inv_stD_cg(chi2,conf,quark->mass,100000,residue,chi1);
	communicate_ev_and_od_color_borders(chi2);
      }
    
    //array to store temp results
    complex *point_result=nissa_malloc("point_result",loc_vol,complex);
    
    /////////////////////// chiral cond and its susceptivity ///////////////////////
    
    //powers of chi: if 1, compute condensate, 2 compute its susceptivity
    for(int ichi=1;ichi<=(comp_susc?2:1);ichi++)
      {
	color **chi=(ichi==1)?chi1:chi2; //if first case use the vector containing the inverse of M, otherwise its square
	
	//summ the prod of EVN and ODD parts
	vector_reset(point_result); //reset the point result
	for(int par=0;par<2;par++) //loop on parity of sites
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh) //loop on sites
	    for(int ic=0;ic<3;ic++) //for every color takes the trace with conjugate of original source
	      complex_summ_the_conj1_prod(point_result[loclx_of_loceo[par][ieo]],rnd[par][ieo][ic],chi[par][ieo][ic]);
	THREAD_BARRIER();
	
	//chir cond: deg/4vol
	complex temp;
	complex_vector_glb_collapse(temp,point_result,loc_vol);
	if(IS_MASTER_THREAD) //normalization: <\bar{\psi}\psi>=Nf/4glb_vol <Tr M^-1>
	  { //susceptivity disconnected: \chi_\psi\psi=-Nf/4glb_vol <Tr M^-2>
	    if(ichi==1) complex_prod_double(putpourri->chiral_cond,temp,quark->deg/(4.0*glb_vol));
	    else        complex_prod_double(putpourri->chiral_cond_susc,temp,-quark->deg/(4.0*glb_vol));
	  }
      }
    
    //compute the derivative of M applied to \chi
    if(comp_susc)
      {
	//generate the source
	for(int par=0;par<2;par++)
	  {
	    vector_reset(app[par]);
            NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	      {
		int idw=loceo_neighdw[par][ieo][0],iup=loceo_neighup[par][ieo][0];
		unsafe_su3_prod_color(app[par][ieo],conf[par][ieo][0],chi1[!par][iup]);
		su3_dag_summ_the_prod_color(app[par][ieo],conf[!par][idw][0],chi1[!par][idw]);
		color_prod_double(app[par][ieo],app[par][ieo],0.5);
	      }
	    set_borders_invalid(app[par]);
	  }	    
	
	//invert
	inv_stD_cg(chi2,conf,quark->mass,100000,residue,app);
	communicate_ev_and_od_color_borders(chi2);
      }
    
    ///////////////////// energy, barionic and pressure density ////////////////
    //compute forward derivative and backward, in turn
    //take into account that backward one must be conjugated
    complex res_fw_bw[4][2];
    for(int mu=0;mu<4;mu++)
      compute_fw_bw_der_mel(res_fw_bw[mu],rnd,conf,mu,chi1,point_result);
    
    //combine forward and backward derivative
    if(IS_MASTER_THREAD)
      {
	//energy density
	complex_subt(putpourri->energy_dens,res_fw_bw[0][0],res_fw_bw[0][1]);
	complex_prodassign_double(putpourri->energy_dens,quark->deg/(4.0*glb_vol)/2);
	//adimensional quark density: <n>=Nf/4vspat <Tr <M^1 dM/d\mu>, 1/2 coming from symm der (mu=pot/T)
	complex_summ(putpourri->quark_dens,res_fw_bw[0][0],res_fw_bw[0][1]);
	complex_prodassign_double(putpourri->quark_dens,quark->deg/(4.0*glb_vol)/2);
	//pressure density
	for(int idir=1;idir<4;idir++)
	  {
	    complex_summassign(putpourri->pressure_dens,res_fw_bw[idir][0]);
	    complex_subtassign(putpourri->pressure_dens,res_fw_bw[idir][1]);
	  }
	complex_prodassign_double(putpourri->pressure_dens,quark->deg/(4.0*glb_vol)/2);
      }
    
    //if needed compute the quark number susceptivity
    if(comp_susc)
      { //adimensional, need to be summed to the energy density! 
	complex res_quark_dens_susc_fw_bw[2];
	compute_fw_bw_der_mel(res_quark_dens_susc_fw_bw,rnd,conf,0,chi2,point_result);
	if(IS_MASTER_THREAD)
	  {
	    complex_summ(putpourri->quark_dens_susc,res_quark_dens_susc_fw_bw[0],res_quark_dens_susc_fw_bw[1]);
	    complex_prodassign_double(putpourri->quark_dens_susc,-quark->deg/(4.0*glb_vol)/2);
	  }
      }
    
    //remove stag phases and u1 field, and automatically barrier before collapsing
    rem_backfield_from_conf(conf,u1b);
    addrem_stagphases_to_eo_conf(conf);
    
    //free automatic synchronizing
    nissa_free(point_result);
    for(int par=0;par<2;par++)
      {
	nissa_free(rnd[par]);
	nissa_free(chi1[par]);
	if(comp_susc)
	  {
	    nissa_free(app[par]);
	    nissa_free(chi2[par]);
	  }
      }
  }
  THREADABLE_FUNCTION_END
  
  //measure the above fermionic putpourri
  void measure_fermionic_putpourri(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(theory_pars.fermionic_putpourri_pars.path,conf_created?"w":"a");
    int comp_susc=theory_pars.fermionic_putpourri_pars.compute_susceptivities;
    
    //measure the putpourri for each quark
    int ncopies=theory_pars.fermionic_putpourri_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
	master_fprintf(file,"%d",iconf);    
	for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
	  {
	    fermionic_putpourri_t putpourri;
	    
	    //loop over hits
	    int nhits=theory_pars.fermionic_putpourri_pars.nhits;
	    for(int hit=0;hit<nhits;hit++)
	      {
		verbosity_lv2_master_printf("Evaluating fermionic putpourri for flavor %d/%d, ncopy %d/%d, nhits %d/%d\n",
					    iflav+1,theory_pars.nflavs,icopy+1,ncopies,hit+1,nhits);
		
		//compute and summ
		fermionic_putpourri_t temp;
		fermionic_putpourri(&temp,conf,theory_pars.backfield[iflav],theory_pars.quark_content+iflav,
				    theory_pars.fermionic_putpourri_pars.residue,
				    comp_susc);
		complex_summassign(putpourri.chiral_cond,temp.chiral_cond);
		if(comp_susc) complex_summassign(putpourri.chiral_cond_susc,temp.chiral_cond_susc);
		complex_summassign(putpourri.energy_dens,temp.energy_dens);
		complex_summassign(putpourri.quark_dens,temp.quark_dens);
		if(comp_susc) complex_summassign(putpourri.quark_dens_susc,temp.quark_dens_susc);
		complex_summassign(putpourri.pressure_dens,temp.pressure_dens);
	      }
	    
	   //write results
	   master_fprintf(file,"\t\t%+016.16lg\t%+016.16lg",putpourri.chiral_cond[RE]/nhits,putpourri.chiral_cond[IM]/nhits);
	   if(comp_susc) master_fprintf(file,"\t%+016.16lg\t%+016.16lg",putpourri.chiral_cond_susc[RE]/nhits,
				 putpourri.chiral_cond_susc[IM]/nhits);
	   master_fprintf(file,"\t%+016.16lg\t%+016.16lg",putpourri.energy_dens[RE]/nhits,putpourri.energy_dens[IM]/nhits);
	   master_fprintf(file,"\t%+016.16lg\t%+016.16lg",putpourri.quark_dens[RE]/nhits,putpourri.quark_dens[IM]/nhits);
	   if(comp_susc) master_fprintf(file,"\t%+016.16lg\t%+016.16lg",putpourri.quark_dens_susc[RE]/nhits,
				 putpourri.quark_dens_susc[IM]/nhits);
	   master_fprintf(file,"\t%+016.16lg\t%+016.16lg",putpourri.pressure_dens[RE]/nhits,
			  putpourri.pressure_dens[IM]/nhits);
	  }
	
	master_fprintf(file,"\n");
      }
    
    //close the file
    if(rank==0) fclose(file);
  }
  
  //compute the magnetization
  THREADABLE_FUNCTION_6ARG(magnetization, complex*,magn, quad_su3**,conf, quad_u1**,u1b, quark_content_t*,quark, double,residue, color**,rnd)
  {
    GET_THREAD_ID();
    
    //fixed to Z magnetization
    int mu=1,nu=2;
    
    //allocate source and propagator
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    
    //we need to store phases
    coords *arg=nissa_malloc("arg",loc_vol+bord_vol,coords);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol+bord_vol)
      get_args_of_one_over_L2_quantization(arg[ivol],ivol,mu,nu);
    
    //array to store magnetization on single site (actually storing backward contrib at displaced site)
    complex *point_magn=nissa_malloc("app",loc_vol,complex);
    vector_reset(point_magn);
    
    //we add stagphases and backfield externally because we need them for derivative
    addrem_stagphases_to_eo_conf(conf);
    add_backfield_to_conf(conf,u1b);
    
    //invert
    inv_stD_cg(chi,conf,quark->mass,100000,residue,rnd);
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
    if(IS_MASTER_THREAD)
      unsafe_complex_prod_idouble(*magn,temp,-quark->deg*2*M_PI*quark->charge/(4.0*glb_vol*2*glb_size[mu]*glb_size[nu]));
    
    //free
    for(int par=0;par<2;par++) nissa_free(chi[par]);
    nissa_free(point_magn);
    nissa_free(arg);
  }
  THREADABLE_FUNCTION_END
  
  //compute the magnetization
  void magnetization(complex *magn,quad_su3 **conf,quad_u1 **u1b,quark_content_t *quark,double residue)
  {
    //allocate source and generate it
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    generate_fully_undiluted_eo_source(rnd,RND_GAUSS,-1);
    
    //call inner function
    magnetization(magn,conf,u1b,quark,residue,rnd);

    for(int par=0;par<2;par++) nissa_free(rnd[par]);
  }
  
  //measure magnetization
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(theory_pars.magnetization_pars.path,conf_created?"w":"a");
    
    int ncopies=theory_pars.magnetization_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
	master_fprintf(file,"%d",iconf);
	
	//measure magnetization for each quark
	for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
	  {
	    complex magn={0,0};
	    
	    //loop over hits
	    int nhits=theory_pars.magnetization_pars.nhits;
	    for(int hit=0;hit<nhits;hit++)
	      {
		verbosity_lv2_master_printf("Evaluating magnetization for flavor %d/%d, ncopies %d/%d nhits %d/%d\n",
					    iflav+1,theory_pars.nflavs,icopy+1,ncopies,hit+1,nhits);
	    
		//compute and summ
		complex temp;
		magnetization(&temp,conf,theory_pars.backfield[iflav],theory_pars.quark_content+iflav,
			      theory_pars.magnetization_pars.residue);
		complex_summ_the_prod_double(magn,temp,1.0/nhits);
	      }
	    
	    //output
	    master_fprintf(file,"\t%+016.16lg \t%+016.16lg",magn[RE],magn[IM]);
	  }
	
	master_fprintf(file,"\n");
      }
    
    if(rank==0) fclose(file);
  }

  //compute the local pseudoscalar correlator in "time" direction (that can be all but time)
  void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created,int dir)
  {
    char dir_name[5]="txyz";
    FILE *file=open_file(theory_pars.pseudo_corr_pars.path,conf_created?"w":"a");
    
    int nflavs=theory_pars.nflavs;
    
    //allocate source
    color *source[2]={nissa_malloc("source_e",loc_volh+bord_volh,color),
		      nissa_malloc("source_o",loc_volh+bord_volh,color)};
    
    //allocate propagators
    color *prop[nflavs][2];
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int EO=0;EO<2;EO++)
	prop[iflav][EO]=nissa_malloc("prop",loc_volh+bord_volh,color);
    
    //allocate local and global contraction
    complex *loc_contr=nissa_malloc("loc_contr",glb_size[dir]*nflavs*(nflavs+1)/2,complex);
    complex *glb_contr=nissa_malloc("glb_contr",glb_size[dir]*nflavs*(nflavs+1)/2,complex);
    vector_reset(loc_contr);
    
    //loop over the hits
    int nhits=theory_pars.pseudo_corr_pars.nhits;
    for(int hit=0;hit<nhits;hit++)
      {
	verbosity_lv2_master_printf("Evaluating pseudoscalar %c correlator, hit %d/%d\n",dir_name[dir],hit+1,nhits);
	
	//generate the source on an even site
	int twall=(int)rnd_get_unif(&glb_rnd_gen,0,glb_size[dir]/2)*2;
	generate_fully_undiluted_eo_source(source,RND_Z4,twall,dir);
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
    glb_reduce_complex_vect(glb_contr,loc_contr,glb_size[dir]*nflavs*(nflavs+1)/2);
    
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
}
