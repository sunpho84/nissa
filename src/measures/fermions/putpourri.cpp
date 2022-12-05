#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"

#include "putpourri.hpp"

#include "stag.hpp"

namespace nissa
{
  using namespace stag;
  
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
  
  //compute the fermionic putpourri for a single conf and hit
  void fermionic_putpourri(fermionic_putpourri_t* putpourri,rnd_t rnd_type,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,quark_content_t* quark,double residue,int comp_susc)
  {
    crash("#warning toredo");
    
    // THREAD_BARRIER();
    
    // //allocate
    // color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    // color *_chi1[2]={nissa_malloc("chi1_EVN",loc_volh+bord_volh,color),nissa_malloc("chi1_ODD",loc_volh+bord_volh,color)};
    // color *_chi2[2],*_app[2];
    // color **chi1=_chi1,**chi2=_chi2,**app=_app;
    // if(comp_susc)
    //   for(int par=0;par<2;par++)
    // 	{
    // 	  chi2[par]=nissa_malloc("chi2_EO",loc_volh+bord_volh,color);
    // 	  app[par]=nissa_malloc("app_EO",loc_volh+bord_volh,color);
    // 	}
    
    // //generate the source
    // generate_fully_undiluted_eo_source(rnd,rnd_type,-1);
    
    // //we add backfield externally because we need them for derivative
    // add_backfield_with_stagphases_to_conf(conf,u1b);
    
    // //invert
    // inv_stD_cg(chi1,conf,quark->mass,100000,residue,rnd);
    // communicate_ev_and_od_color_borders(chi1);
    
    // //invert for chi2 if susceptivity needed
    // if(comp_susc)
    //   {
    // 	inv_stD_cg(chi2,conf,quark->mass,100000,residue,chi1);
    // 	communicate_ev_and_od_color_borders(chi2);
    //   }
    
    // //array to store temp results
    // complex *point_result=nissa_malloc("point_result",loc_vol,complex);
    
    // /////////////////////// chiral cond and its susceptivity ///////////////////////
    
    // //powers of chi: if 1, compute condensate, 2 compute its susceptivity
    // for(int ichi=1;ichi<=(comp_susc?2:1);ichi++)
    //   {
    // 	color **chi=(ichi==1)?chi1:chi2; //if first case use the vector containing the inverse of M, otherwise its square
	
    // 	//summ the prod of EVN and ODD parts
    // 	vector_reset(point_result); //reset the point result
    // 	for(int par=0;par<2;par++) //loop on parity of sites
    // 	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh) //loop on sites
    // 	    for(int ic=0;ic<NCOL;ic++) //for every color takes the trace with conjugate of original source
    // 	      complex_summ_the_conj1_prod(point_result[loclx_of_loceo[par][ieo]],rnd[par][ieo][ic],chi[par][ieo][ic]);
    // 	  NISSA_PARALLEL_LOOP_END;
    // 	THREAD_BARRIER();
	
    // 	//chir cond: deg/4vol
    // 	complex temp;
    // 	complex_vector_glb_collapse(temp,point_result,loc_vol);
    // 	//DEB_STAG("ichi=%d temp=(%lg %lg)\n",ichi,temp[RE],temp[IM]);
	
    // 	if(IS_MASTER_THREAD) //normalization: <\bar{\psi}\psi>=Nf/4glb_vol <Tr M^-1>
    // 	  { //susceptivity disconnected: \chi_\psi\psi=-Nf/4glb_vol <Tr M^-2>
    // 	    if(ichi==1) complex_prod_double(putpourri->chiral_cond,temp,quark->deg/(4.0*glb_vol));
    // 	    else        complex_prod_double(putpourri->chiral_cond_susc,temp,-quark->deg/(4.0*glb_vol));
    // 	  }
    // 	THREAD_BARRIER();
    //   }
    
    // //compute the derivative of M applied to \chi
    // if(comp_susc)
    //   {
    // 	//generate the source
    // 	for(int par=0;par<2;par++)
    // 	  {
    // 	    vector_reset(app[par]);
    //         NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
    // 	      {
    // 		int idw=loceo_neighdw[par][ieo][0],iup=loceo_neighup[par][ieo][0];
    // 		unsafe_su3_prod_color(app[par][ieo],conf[par][ieo][0],chi1[!par][iup]);
    // 		su3_dag_summ_the_prod_color(app[par][ieo],conf[!par][idw][0],chi1[!par][idw]);
    // 		color_prod_double(app[par][ieo],app[par][ieo],0.5);
    // 	      }
    // 	    NISSA_PARALLEL_LOOP_END;
    // 	    set_borders_invalid(app[par]);
    // 	  }
	
    // 	//invert
    // 	inv_stD_cg(chi2,conf,quark->mass,100000,residue,app);
    // 	communicate_ev_and_od_color_borders(chi2);
    //   }
    
    // ///////////////////// energy, barionic and pressure density ////////////////
    // //compute forward derivative and backward, in turn
    // //take into account that backward one must be conjugated
    // complex res_fw_bw[NDIM][2];
    // for(int mu=0;mu<4;mu++)
    //   compute_fw_bw_der_mel(res_fw_bw[mu],rnd,conf,mu,chi1,point_result);
    
    // //combine forward and backward derivative
    // if(IS_MASTER_THREAD)
    //   {
    // 	//energy density
    // 	complex_subt(putpourri->energy_dens,res_fw_bw[0][0],res_fw_bw[0][1]);
    // 	complex_prodassign_double(putpourri->energy_dens,quark->deg/(4.0*glb_vol)/2);
    // 	//adimensional quark density: <n>=Nf/4vspat <Tr <M^1 dM/d\mu>, 1/2 coming from symm der (mu=pot/T)
    // 	complex_summ(putpourri->quark_dens,res_fw_bw[0][0],res_fw_bw[0][1]);
    // 	complex_prodassign_double(putpourri->quark_dens,quark->deg/(4.0*glb_vol)/2);
    // 	//pressure density
    // 	for(int idir=1;idir<4;idir++)
    // 	  {
    // 	    complex_summassign(putpourri->pressure_dens,res_fw_bw[idir][0]);
    // 	    complex_subtassign(putpourri->pressure_dens,res_fw_bw[idir][1]);
    // 	  }
    // 	complex_prodassign_double(putpourri->pressure_dens,quark->deg/(4.0*glb_vol)/2);
    //   }
    // THREAD_BARRIER();
    
    // //if needed compute the quark number susceptibility
    // if(comp_susc)
    //   { //adimensional, need to be summed to the energy density!
    // 	complex res_quark_dens_susc_fw_bw[2];
    // 	compute_fw_bw_der_mel(res_quark_dens_susc_fw_bw,rnd,conf,0,chi2,point_result);
    // 	if(IS_MASTER_THREAD)
    // 	  {
    // 	    complex_summ(putpourri->quark_dens_susc,res_quark_dens_susc_fw_bw[0],res_quark_dens_susc_fw_bw[1]);
    // 	    complex_prodassign_double(putpourri->quark_dens_susc,-quark->deg/(4.0*glb_vol)/2);
    // 	  }
    // 	THREAD_BARRIER();
    //   }
    
    // //remove stag phases and u1 field, and automatically barrier before collapsing
    // rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    // //free automatic synchronizing
    // nissa_free(point_result);
    // for(int par=0;par<2;par++)
    //   {
    // 	nissa_free(rnd[par]);
    // 	nissa_free(chi1[par]);
    // 	if(comp_susc)
    // 	  {
    // 	    nissa_free(app[par]);
    // 	    nissa_free(chi2[par]);
    // 	  }
    //   }
  }
  
  //measure the above fermionic putpourri
  void measure_fermionic_putpourri(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,fermionic_putpourri_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    crash("reimplement");
    
    // FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    // int comp_susc=meas_pars.compute_susc;
    
    // //measure the putpourri for each quark
    // int ncopies=meas_pars.ncopies;
    // for(int icopy=0;icopy<ncopies;icopy++)
    //   {
    // 	master_fprintf(file,"%d",iconf);
    // 	for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
    // 	  {
    // 	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
    // 	    fermionic_putpourri_t putpourri;
	    
    // 	    //loop over hits
    // 	    int nhits=meas_pars.nhits;
    // 	    for(int hit=0;hit<nhits;hit++)
    // 	      {
    // 		verbosity_lv2_master_printf("Evaluating fermionic putpourri for flavor %d/%d, ncopy %d/%d, nhits %d/%d\n",
    // 					    iflav+1,theory_pars.nflavs(),icopy+1,ncopies,hit+1,nhits);
		
    // 		//compute and summ
    // 		fermionic_putpourri_t temp;
    // 		fermionic_putpourri(&temp,meas_pars.rnd_type,conf,theory_pars.backfield[iflav],&theory_pars.quarks[iflav],meas_pars.residue,comp_susc);
    // 		complex_summassign(putpourri.chiral_cond,temp.chiral_cond);
    // 		if(comp_susc) complex_summassign(putpourri.chiral_cond_susc,temp.chiral_cond_susc);
    // 		complex_summassign(putpourri.energy_dens,temp.energy_dens);
    // 		complex_summassign(putpourri.quark_dens,temp.quark_dens);
    // 		if(comp_susc) complex_summassign(putpourri.quark_dens_susc,temp.quark_dens_susc);
    // 		complex_summassign(putpourri.pressure_dens,temp.pressure_dens);
    // 	      }
	    
    // 	   //write results
    // 	   master_fprintf(file,"\t\t%+16.16lg\t%+16.16lg",putpourri.chiral_cond[RE]/nhits,putpourri.chiral_cond[IM]/nhits);
    // 	   if(comp_susc) master_fprintf(file,"\t%+16.16lg\t%+16.16lg",putpourri.chiral_cond_susc[RE]/nhits,
    // 				 putpourri.chiral_cond_susc[IM]/nhits);
    // 	   master_fprintf(file,"\t%+16.16lg\t%+16.16lg",putpourri.energy_dens[RE]/nhits,putpourri.energy_dens[IM]/nhits);
    // 	   master_fprintf(file,"\t%+16.16lg\t%+16.16lg",putpourri.quark_dens[RE]/nhits,putpourri.quark_dens[IM]/nhits);
    // 	   if(comp_susc) master_fprintf(file,"\t%+16.16lg\t%+16.16lg",putpourri.quark_dens_susc[RE]/nhits,
    // 				 putpourri.quark_dens_susc[IM]/nhits);
    // 	   master_fprintf(file,"\t%+16.16lg\t%+16.16lg",putpourri.pressure_dens[RE]/nhits,
    // 			  putpourri.pressure_dens[IM]/nhits);
    // 	  }
	
    // 	master_fprintf(file,"\n");
    //   }
    
    // //close the file
    // if(rank==0) fclose(file);
  }
  
  //fermionic putpourri
  std::string fermionic_putpourri_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasPutpourri\n";
    if(is_nonstandard()||full) os<<base_fermionic_meas_t::get_str(full);
    if(compute_susc!=def_compute_susc()||full)  os<<" ComputeSusc\t=\t"<<compute_susc<<"\n";
    
    return os.str();
  }
}
