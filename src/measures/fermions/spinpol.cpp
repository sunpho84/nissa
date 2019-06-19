#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "hmc/theory_pars.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "measures/fermions/stag.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "measures/gauge/topological_charge.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "operations/smearing/recursive_Wflower.hpp"
#include "spinpol.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  using namespace stag;
  typedef std::vector<std::pair<int,int> > op_list_t;
  
  namespace
  {
    const int nPHIETA=2,PHI=0,ETA=1;
    int ncopies,nflavs,nhits,nmeas,nops;
    int ind_copy_hit(int icopy,int ihit){return ihit+nhits*icopy;}
    int ind_copy_flav_hit_meas(int icopy,int iflav,int ihit,int imeas){return imeas+nmeas*(ihit+nhits*(iflav+nflavs*icopy));}
    int ind_copy_flav_hit_phieta(int icopy,int iflav,int ihit,int iPHIETA){return iPHIETA+nPHIETA*(ihit+nhits*(iflav+nflavs*icopy));}
  }
  
  //make the complex-double product
  void compute_tens_dens_topo_correlation(complex *spinpol_dens,complex *tens_dens,double *topo_dens)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) complex_prod_double(spinpol_dens[ivol],tens_dens[ivol],topo_dens[ivol]);
    THREAD_BARRIER();
  }
  
  //compute the spin-polarization for all flavors
  THREADABLE_FUNCTION_5ARG(measure_spinpol, theory_pars_t*,tp, spinpol_meas_pars_t*,mp,int,iconf, int,conf_created, quad_su3**,glu_conf)
  {
    verbosity_lv1_master_printf("Evaluating spinpol\n");
    
    //set-up the smoother
    smooth_pars_t &sp=mp->smooth_pars;
    if(sp.method!=smooth_pars_t::WFLOW) crash("spinpol makes sense only with Wflow");
    Wflow_pars_t &Wf=sp.Wflow;
    int nflows=Wf.nflows;
    double dt=Wf.dt;
    int meas_each=sp.meas_each_nsmooth;
    
    //take number of flavors and operators
    nflavs=tp->nflavs();
    nops=mp->nops();
    
    //open the file
    FILE *fout=open_file(mp->path,conf_created?"w":"a");
    
    //allocate tens and spinpol
    complex *tens_dens=nissa_malloc("tens_dens",loc_vol+bord_vol,complex);
    complex *spinpol_dens=nissa_malloc("spinpol_dens",loc_vol,complex);
    //allocate point and local results
    double *topo_dens[2];
    for(int i=0;i<2;i++) topo_dens[i]=nissa_malloc("topo_dens",loc_vol,double);
    //operator applied to a field
    color *chiop[2]={nissa_malloc("chiop_EVN",loc_volh+bord_volh,color),nissa_malloc("chiop_ODD",loc_volh+bord_volh,color)};
    //temporary vectors
    color *temp[2][2];
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	temp[itemp][eo]=nissa_malloc("temp",loc_volh+bord_volh,color);
    
    //Create the mask
    int mask[nops],shift[nops];
    for(int iop=0;iop<nops;iop++)
      {
	int spin=mp->operators[iop].first;
	int taste=mp->operators[iop].second;
	shift[iop]=(spin^taste);
	mask[iop]=form_stag_op_pattern(spin,taste);
	verbosity_lv3_master_printf(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
      }
    
    //allocate the smoothed conf
    quad_su3 *smoothed_conf=nissa_malloc("smoothed_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_vector(smoothed_conf,glu_conf);
    //allocate the fermion (possibly stouted) conf
    quad_su3 *ferm_conf[2];
    for(int eo=0;eo<2;eo++) ferm_conf[eo]=nissa_malloc("ferm_conf",loc_volh+bord_volh+edge_volh,quad_su3);
    
    //allocate sources, to be flown
    ncopies=mp->ncopies;
    nmeas=mp->smooth_pars.nmeas_nonzero()+1;
    nhits=mp->nhits;
    
    if(mp->use_adjoint_flow)
      {
	color *temp_eta[2]={nissa_malloc("eta_EVN",loc_volh+bord_volh,color),nissa_malloc("eta_ODD",loc_volh+bord_volh,color)};
	color *temp_phi[2]={nissa_malloc("phi_EVN",loc_volh+bord_volh,color),nissa_malloc("phi_ODD",loc_volh+bord_volh,color)};
	
	int ntot_eta=ind_copy_hit(ncopies-1,nhits-1)+1;
	int ntot_phi=ind_copy_flav_hit_meas(ncopies-1,nflavs-1,nhits-1,nmeas-1)+1;
	color *eta[ntot_eta];
	color *phi[ntot_phi];
	for(int is=0;is<ntot_eta;is++) eta[is]=nissa_malloc("eta",loc_vol+bord_vol,color);
	for(int is=0;is<ntot_phi;is++) phi[is]=nissa_malloc("phi",loc_vol+bord_vol,color);
	
	//fill all the sources
	for(int icopy=0;icopy<ncopies;icopy++)
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      int ieta=ind_copy_hit(icopy,ihit);
	      generate_fully_undiluted_lx_source(eta[ieta],mp->rnd_type,-1);
	      
	      for(int iflav=0;iflav<nflavs;iflav++)
		for(int imeas=0;imeas<nmeas;imeas++)
		  vector_copy(phi[ind_copy_flav_hit_meas(icopy,iflav,ihit,imeas)],eta[ieta]);
	    }
	
	//the reecursive flower, need to cache backward integration
	paste_eo_parts_into_lx_vector(smoothed_conf,glu_conf);
	recursive_Wflower_t recu(Wf,smoothed_conf);
	//the adjoint flower needed for fermionic source
	fermion_adjoint_flower_t<> adj_ferm_flower(dt,all_dirs);
	
	// int nevol=0;
	// for(int i=sp.nsmooth();i>=0;i--)
	//   {
	// 	master_printf("\n");
	// 	nevol+=recu.update(i);
	//   }
	// master_printf("nevol: %d\n",nevol);
	
	//compute the topological charge and the product of topological and tensorial density
	
	//1) draw all hits and fill them at all time we want to measure, evolve to time 0
	//2) use plain guon action from t 0 (eq.7.7 of 1302.5246) to evolve back from tmeas to 0
	//3) solve all hits using fermion action and plug all operators
	//4) take the trace with the hits
	//5) compute topocharge at all intermediate times
	
	//at each step it goes from iflow+1 to iflow
	for(int iflow=nflows-1;iflow>=0;iflow--)
	  {
	    //update conf to iflow
	    double t=dt*iflow;
	    verbosity_lv2_master_printf(" flow back to %d/%d, t %lg\n",iflow,nflows,t);
	    recu.update(iflow);
	    
	    //make the flower generate the intermediate step between iflow and iflow+1
	    adj_ferm_flower.generate_intermediate_steps(smoothed_conf);
	    
	    //have to flow back all sources for which iflow is smaller than meas_each*imeas
	    int imeas_min=iflow/meas_each+1;
	    for(int iflav=0;iflav<nflavs;iflav++)
	      {
		adj_ferm_flower.add_or_rem_backfield_to_confs(0,tp->backfield[iflav]);
		for(int imeas=imeas_min;imeas<nmeas;imeas++)
		  for(int icopy=0;icopy<ncopies;icopy++)
		    for(int ihit=0;ihit<nhits;ihit++)
		      adj_ferm_flower.flow_fermion(phi[ind_copy_flav_hit_meas(icopy,iflav,ihit,imeas)]);
		adj_ferm_flower.add_or_rem_backfield_to_confs(1,tp->backfield[iflav]);
	      }
	  }
	
	//build fermionic conf from gluonic one
	for(int eo=0;eo<2;eo++) vector_copy(ferm_conf[eo],glu_conf[eo]);
	if(tp->stout_pars.nlevels)
	  stout_smear(ferm_conf,ferm_conf,&tp->stout_pars);
	
	//put the operator
	for(int imeas=0;imeas<nmeas;imeas++)
	  for(int icopy=0;icopy<ncopies;icopy++)
	    for(int ihit=0;ihit<nhits;ihit++)
	      for(int iflav=0;iflav<nflavs;iflav++)
		{
		  split_lx_vector_into_eo_parts(temp_eta,phi[ind_copy_flav_hit_meas(icopy,iflav,ihit,imeas)]);
		  mult_Minv(temp_phi,ferm_conf,tp,iflav,mp->residue,temp_eta);
		  paste_eo_parts_into_lx_vector(phi[ind_copy_flav_hit_meas(icopy,iflav,ihit,imeas)],temp_phi);
		}
      	
	fermion_flower_t<color,4> ferm_flower(dt,all_dirs);
	for(int iflow=0;iflow<=nflows;iflow++)
	  {
	    //take current meas index
	    int imeas=iflow/meas_each;
	    
	    //measure only if flow is a measure time
	    if(imeas*meas_each==iflow)
	      {
		//create the fermionic conf
		split_lx_vector_into_eo_parts(ferm_conf,smoothed_conf);
		if(tp->stout_pars.nlevels)
		  stout_smear(ferm_conf,ferm_conf,&tp->stout_pars);
		
		//plaquette and local charge
		double plaq[2];
		double tot_charge[2];
		double tot_charge2[2];
		
		//use gauge or ferm_conf in alternartive
		for(int gauge_ferm_conf=0;gauge_ferm_conf<2;gauge_ferm_conf++)
		  {
		    quad_su3 *topo_conf;
		    if(gauge_ferm_conf==0) topo_conf=smoothed_conf;
		    else
		      {
			topo_conf=nissa_malloc("TempConf",loc_vol+bord_vol+edge_vol,quad_su3);
			paste_eo_parts_into_lx_vector(topo_conf,ferm_conf);
		      }
		    plaq[gauge_ferm_conf]=global_plaquette_lx_conf(topo_conf);
		    local_topological_charge(topo_dens[gauge_ferm_conf],topo_conf);
		    if(gauge_ferm_conf) nissa_free(topo_conf);
		    
		    //total topological charge
		    double_vector_glb_collapse(&tot_charge[gauge_ferm_conf],topo_dens[gauge_ferm_conf],loc_vol);
		    tot_charge2[gauge_ferm_conf]=double_vector_glb_norm2(topo_dens[gauge_ferm_conf],loc_vol);
		  }
		
		for(int icopy=0;icopy<ncopies;icopy++)
		  for(int iflav=0;iflav<nflavs;iflav++)
		    for(int iop=0;iop<nops;iop++)
		      {
			//compute the local tensorial density
			vector_reset(tens_dens);
			for(int ihit=0;ihit<nhits;ihit++)
			  {
			    split_lx_vector_into_eo_parts(temp_phi,phi[ind_copy_flav_hit_meas(icopy,iflav,ihit,imeas)]);
			    split_lx_vector_into_eo_parts(temp_eta,eta[ind_copy_hit(icopy,ihit)]);
			    
			    summ_dens(tens_dens,chiop,temp[0],temp[1],ferm_conf,tp->backfield[iflav],shift[iop],mask[iop],temp_phi,temp_eta);
			  }
			
			//compute the average tensorial density
			complex tens;
			double_vector_prodassign_double((double*)tens_dens,1.0/(glb_vol*nhits),loc_vol*2);
			complex_vector_glb_collapse(tens,tens_dens,loc_vol);
			
			//compute correlation with topocharge
			complex spinpol;
			for(int gauge_ferm_conf=0;gauge_ferm_conf<2;gauge_ferm_conf++)
			  {
			    compute_tens_dens_topo_correlation(spinpol_dens,tens_dens,topo_dens[gauge_ferm_conf]);
			    complex_vector_glb_collapse(spinpol,spinpol_dens,loc_vol);
			    
			    master_fprintf(fout,"%d\t",iconf);
			    master_fprintf(fout,"%d\t",icopy);
			    master_fprintf(fout,"%d\t",imeas*meas_each);
			    master_fprintf(fout,"%d\t",iflav);
			    master_fprintf(fout,"%d,%d\t",mp->operators[iop].first,mp->operators[iop].second);
			    master_fprintf(fout,"%d\t",gauge_ferm_conf);
			    master_fprintf(fout,"%+16.16lg\t",plaq[gauge_ferm_conf]);
			    master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",tot_charge[gauge_ferm_conf],tot_charge2[gauge_ferm_conf]);
			    master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",spinpol[RE],spinpol[IM]);
			    master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",tens[RE],tens[IM]);
			    master_fprintf(fout,"\n");
			  }
		      }
	      }
	    
	    //update conf to iflow
	    double t=dt*iflow;
	    //verbosity_lv2_
	    master_printf(" flow forward to %d/%d, t %lg, initial plaquette: %.16lg\n",iflow,nflows,t,global_plaquette_lx_conf(smoothed_conf));
	    
	    //make the flower generate the intermediate step between iflow-1 and iflow
	    ferm_flower.generate_intermediate_steps(smoothed_conf);
	    for(int iflav=0;iflav<nflavs;iflav++)
	      {
	    	ferm_flower.add_or_rem_backfield_to_confs(0,tp->backfield[iflav]);
	    	for(int icopy=0;icopy<ncopies;icopy++)
	    	  for(int ihit=0;ihit<nhits;ihit++)
		    for(int jmeas=imeas;jmeas<nmeas;jmeas++)
		      ferm_flower.flow_fermion(phi[ind_copy_flav_hit_meas(icopy,iflav,ihit,jmeas)]);
		ferm_flower.add_or_rem_backfield_to_confs(1,tp->backfield[iflav]);
	      }
	    ferm_flower.prepare_for_next_flow(smoothed_conf);
	  }
	
	//free
	for(int i=0;i<ntot_phi;i++) nissa_free(phi[i]);
	for(int i=0;i<ntot_eta;i++) nissa_free(eta[i]);
	for(int eo=0;eo<2;eo++)
	  {
	    nissa_free(temp_eta[eo]);
	    nissa_free(temp_phi[eo]);
	  }
      }
    else
      {
	//Allocate and create all fields, sources and prop*source.
	//fields must be accessed though the index, the last component
	//decides wheter the field is a source (when using value
	//"ieta") or the result
	int nfields=ind_copy_flav_hit_phieta(ncopies-1,nflavs-1,nhits-1,nPHIETA-1)+1;
	color *fields[nfields][2];
	for(int ifield=0;ifield<nfields;ifield++)
	  for(int eo=0;eo<2;eo++)
	    fields[ifield][eo]=nissa_malloc("field",loc_volh+bord_volh,color);
	
	//create the fermionic conf
	for(int eo=0;eo<2;eo++) vector_copy(ferm_conf[eo],glu_conf[eo]);
	if(tp->stout_pars.nlevels)
	  stout_smear(ferm_conf,ferm_conf,&tp->stout_pars);
	
	color *temp_flow=nissa_malloc("temp_flow",loc_vol+bord_vol,color);
	
	for(int icopy=0;icopy<ncopies;icopy++)
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      int isource=ind_copy_flav_hit_phieta(icopy,0,ihit,ETA);
	      generate_fully_undiluted_eo_source(fields[isource],mp->rnd_type,-1);
	      
	      for(int iflav=0;iflav<nflavs;iflav++)
		{
		  int ieta=ind_copy_flav_hit_phieta(icopy,iflav,ihit,ETA);
		  int iphi=ind_copy_flav_hit_phieta(icopy,iflav,ihit,PHI);
		  
		  //if not first flavour, copy the source, and split it
		  if(iflav!=0)
		    for(int eo=0;eo<2;eo++)
		      vector_copy(fields[ieta][eo],fields[isource][eo]);
		  mult_Minv(fields[iphi],ferm_conf,tp,iflav,mp->residue,fields[ieta]);
		}
	    }
	
	fermion_flower_t<color,4> ferm_flower(dt,all_dirs);
	for(int iflow=0;iflow<=nflows;iflow++)
	  {
	    //take current meas index
	    int imeas=iflow/meas_each;
	    
	    if(imeas*meas_each==iflow)
	      {
		//create the fermionic conf
		split_lx_vector_into_eo_parts(ferm_conf,smoothed_conf);
		if(tp->stout_pars.nlevels)
		  stout_smear(ferm_conf,ferm_conf,&tp->stout_pars);
		
		//plaquette and local charge
		//plaquette and local charge
		double plaq[2];
		double tot_charge[2];
		double tot_charge2[2];
		
		//use gauge or ferm_conf in alternartive
		for(int gauge_ferm_conf=0;gauge_ferm_conf<2;gauge_ferm_conf++)
		  {
		    quad_su3 *topo_conf;
		    if(gauge_ferm_conf==0) topo_conf=smoothed_conf;
		    else
		      {
			topo_conf=nissa_malloc("TempConf",loc_vol+bord_vol+edge_vol,quad_su3);
			paste_eo_parts_into_lx_vector(topo_conf,ferm_conf);
		      }
		    plaq[gauge_ferm_conf]=global_plaquette_lx_conf(topo_conf);
		    local_topological_charge(topo_dens[gauge_ferm_conf],topo_conf);
		    if(gauge_ferm_conf) nissa_free(topo_conf);
		    
		    //total topological charge
		    double_vector_glb_collapse(&tot_charge[gauge_ferm_conf],topo_dens[gauge_ferm_conf],loc_vol);
		    tot_charge2[gauge_ferm_conf]=double_vector_glb_norm2(topo_dens[gauge_ferm_conf],loc_vol);
		  }
		
		for(int icopy=0;icopy<ncopies;icopy++)
		  for(int iflav=0;iflav<nflavs;iflav++)
		    for(int iop=0;iop<nops;iop++)
		      {
			//compute the local tensorial density
			vector_reset(tens_dens);
			for(int ihit=0;ihit<nhits;ihit++)
			  {
			    int ieta=ind_copy_flav_hit_phieta(icopy,iflav,ihit,ETA);
			    int iphi=ind_copy_flav_hit_phieta(icopy,iflav,ihit,PHI);
			    
			    summ_dens(tens_dens,chiop,temp[0],temp[1],ferm_conf,tp->backfield[iflav],shift[iop],mask[iop],fields[iphi],fields[ieta]);
			  }
			
			//compute the average tensorial density
			complex tens;
			double_vector_prodassign_double((double*)tens_dens,1.0/(glb_vol*nhits),loc_vol*2);
			complex_vector_glb_collapse(tens,tens_dens,loc_vol);
			
			//compute correlation with topocharge
			complex spinpol;
			for(int gauge_ferm_conf=0;gauge_ferm_conf<2;gauge_ferm_conf++)
			  {
			    compute_tens_dens_topo_correlation(spinpol_dens,tens_dens,topo_dens[gauge_ferm_conf]);
			    complex_vector_glb_collapse(spinpol,spinpol_dens,loc_vol);
			    
			    master_fprintf(fout,"%d\t",iconf);
			    master_fprintf(fout,"%d\t",icopy);
			    master_fprintf(fout,"%d\t",imeas*meas_each);
			    master_fprintf(fout,"%d\t",iflav);
			    master_fprintf(fout,"%d,%d\t",mp->operators[iop].first,mp->operators[iop].second);
			    master_fprintf(fout,"%d\t",gauge_ferm_conf);
			    master_fprintf(fout,"%+16.16lg\t",plaq[gauge_ferm_conf]);
			    master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",tot_charge[gauge_ferm_conf],tot_charge2[gauge_ferm_conf]);
			    master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",spinpol[RE],spinpol[IM]);
			    master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",tens[RE],tens[IM]);
			    master_fprintf(fout,"\n");
			  }
		      }
	      }
	    
	    //update conf to iflow
	    double t=dt*iflow;
	    //verbosity_lv2_
	    master_printf(" flow forward to %d/%d, t %lg, initial plaquette: %.16lg\n",iflow,nflows,t,global_plaquette_lx_conf(smoothed_conf));
	    
	    //make the flower generate the intermediate step between iflow-1 and iflow
	    ferm_flower.generate_intermediate_steps(smoothed_conf);
	    for(int iflav=0;iflav<nflavs;iflav++)
	      {
	    	ferm_flower.add_or_rem_backfield_to_confs(0,tp->backfield[iflav]);
	    	for(int icopy=0;icopy<ncopies;icopy++)
	    	  for(int ihit=0;ihit<nhits;ihit++)
	    	    for(int iPHIETA=0;iPHIETA<nPHIETA;iPHIETA++)
		      {
			paste_eo_parts_into_lx_vector(temp_flow,fields[ind_copy_flav_hit_phieta(icopy,iflav,ihit,iPHIETA)]);
			ferm_flower.flow_fermion(temp_flow);
			split_lx_vector_into_eo_parts(fields[ind_copy_flav_hit_phieta(icopy,iflav,ihit,iPHIETA)],temp_flow);
		      }
	    	ferm_flower.add_or_rem_backfield_to_confs(1,tp->backfield[iflav]);
	      }
	    ferm_flower.prepare_for_next_flow(smoothed_conf);
	  }
	
	for(int ifield=0;ifield<nfields;ifield++)
	  for(int eo=0;eo<2;eo++)
	    nissa_free(fields[ifield][eo]);
	nissa_free(temp_flow);
      }
    
    for(int eo=0;eo<2;eo++)
      {
	nissa_free(chiop[eo]);
	for(int itemp=0;itemp<2;itemp++)
	  nissa_free(temp[itemp][eo]);
	nissa_free(ferm_conf[eo]);
      }
    nissa_free(spinpol_dens);
    nissa_free(tens_dens);
    nissa_free(smoothed_conf);
    for(int gauge_ferm_conf=0;gauge_ferm_conf<2;gauge_ferm_conf++)
      nissa_free(topo_dens[gauge_ferm_conf]);
    
    //close
    close_file(fout);
  }
  THREADABLE_FUNCTION_END
  
  std::string spinpol_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasSpinPol\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(use_adjoint_flow!=def_use_adjoint_flow() or full) os<<" UseAdjointFlow\t=\t"<<use_adjoint_flow<<"\n";
    if(use_ferm_conf_for_gluons!=def_use_ferm_conf_for_gluons() or full) os<<" UseFermConfForGluons\t=\t"<<use_ferm_conf_for_gluons<<"\n";
    if(operators.size())
      {
	os<<" Operators\t=\t{";
	for(size_t i=0;i<operators.size();i++)
	  {
	    os<<"("<<operators[i].first<<","<<operators[i].second<<")";
	    if(i!=operators.size()-1) os<<",";
	  }
	os<<"}\n";
      }
    os<<smooth_pars.get_str(full);
    
    return os.str();
  }
}
