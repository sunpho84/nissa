#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "hmc/theory_pars.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "operations/smearing/recursive_Wflower.hpp"
#include "spinpol.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  using namespace stag;	//In namespace stag there are all functions of mesons.cpp (L)
  typedef std::vector<std::pair<int,int> > op_list_t;
  
  namespace
  {
    int nhits,nops;
    int ind_meas_hit(int imeas,int ihit){return ihit+nhits*imeas;}
    int ind_op_flav(int iop,int iflav){return iop+nops*iflav;}
  }
  
  //compute the tensorial density
  void summ_tens_dens(complex *spinpol_dens,color *quark[2],color *temp[2][2],quad_su3 *ferm_conf[2],quad_u1 *backfield[2],int shift,int mask,color *eta[2])
  {
    GET_THREAD_ID();
    
    apply_op(quark,temp[0],temp[1],ferm_conf,backfield,shift,eta);
    put_stag_phases(quark,mask);
    
    for(int eo=0;eo<2;eo++)
      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	{
	  int ivol=loclx_of_loceo[eo][ieo];
	  complex prod;
	  color_scalar_prod(prod,eta[eo][ieo],quark[eo][ieo]);
	  complex_summassign(spinpol_dens[ivol],prod);
	}
    THREAD_BARRIER();
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
    int nflavs=tp->nflavs();
    nops=mp->nops();
    
    //open the file
    FILE *fout=open_file(mp->path,conf_created?"w":"a");
    
    //allocate point and local results
    double *topo_dens=nissa_malloc("topo_dens",loc_vol,double);
    complex *spinpol_dens[nflavs*nops],spinpol[nflavs*nops];
    for(int iflav_op=0;iflav_op<nflavs*nops;iflav_op++) spinpol_dens[iflav_op]=nissa_malloc("spinpol_dens",loc_vol,complex);
    complex tens[nflavs*nops];
    //tens density
    complex *tens_dens[nflavs*nops];
    for(int iflav_op=0;iflav_op<nflavs*nops;iflav_op++) tens_dens[iflav_op]=nissa_malloc("tens_dens",loc_vol+bord_vol,complex);
    
    //allocate Dirac equation solution and source
    color *eta[2]={nissa_malloc("eta_EVN",loc_volh+bord_volh,color),nissa_malloc("eta_ODD",loc_volh+bord_volh,color)};
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    color *chiop[2]={nissa_malloc("chiop_EVN",loc_volh+bord_volh,color),nissa_malloc("chiop_ODD",loc_volh+bord_volh,color)};
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
    quad_su3 *smoothed_conf=nissa_malloc("smoothed_conf",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_vector(smoothed_conf,glu_conf);
    //allocate the fermion (possibly stouted) conf
    quad_su3 *ferm_conf[2];
    for(int eo=0;eo<2;eo++) ferm_conf[eo]=nissa_malloc("ferm_conf",loc_volh+bord_volh+edge_volh,quad_su3);
    
    //allocate sources, to be flown
    int nmeas=mp->smooth_pars.nmeas_nonzero()+1;
    nhits=mp->nhits;
    color *source[nmeas*nhits];
    for(int is=0;is<nmeas*nhits;is++) source[is]=nissa_malloc("source",loc_vol+bord_vol,color);
    
    //fill all the sources, putting for all measures the same hit
    for(int ihit=0;ihit<nhits;ihit++)
      {
	color *source0=source[ind_meas_hit(0/*imeas*/,ihit)];
	generate_fully_undiluted_lx_source(source0,RND_Z4,-1);
	for(int imeas=1;imeas<nmeas;imeas++) vector_copy(source[ind_meas_hit(imeas,ihit)],source0);
      }
    
    //the reecursive flower, need to cache backward integration
    recursive_Wflower_t recu(Wf,smoothed_conf);
    //the afjoint flower needed for fermionic source
    fermion_adjoint_flower_t ferm_flower(dt,all_dirs,true);
    
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
    
    int ncopies=mp->ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
	verbosity_lv2_master_printf("Copy: %d/%d\n",icopy,ncopies);
	
	for(int iflow=nflows;iflow>=0;iflow--)
	  {
	    //update conf
	    double t=dt*iflow;
	    verbosity_lv2_master_printf(" iflow %d/%d, t %lg\n",iflow,nflows,t);
	    recu.update(iflow);
	    
	    //make the flower generate the intermediate step
	    ferm_flower.generate_intermediate_steps(smoothed_conf);
	    
	    //have to flow back all sources for which iflow is smaller than meas_each*imeas
	    int imeas_min=nflows/meas_each;
	    for(int imeas=imeas_min;imeas<nmeas;imeas++)
	      for(int ihit=0;ihit<nhits;ihit++)
		ferm_flower.flow_fermion(source[ind_meas_hit(imeas,ihit)]);
	    
	    //check if this is a meas_each
	    int imeas=iflow/meas_each;
	    printf(" iflow: %d meas_each: %d imeas: %d imeas*meas_each==iflow: %d\n",iflow,meas_each,imeas,imeas*meas_each==iflow);
	    if(imeas*meas_each==iflow)
	      {
		//build fermionic conf
		split_lx_vector_into_eo_parts(ferm_conf,smoothed_conf);
		if(tp->stout_pars.nlevels)
		  stout_smear(ferm_conf,ferm_conf,&tp->stout_pars);
		
		//reset the local density of topological charge
		for(int iflav=0;iflav<nflavs;iflav++)
		  for(int iop=0;iop<nops;iop++)
		    vector_reset(tens_dens[ind_op_flav(iop,iflav)]);
		
		//evaluate the tensorial density for all quarks
		for(int ihit=0;ihit<nhits;ihit++)
		  {
		    int iso=ind_meas_hit(imeas,ihit);
		    split_lx_vector_into_eo_parts(eta,source[iso]);
		    for(int iflav=0;iflav<nflavs;iflav++)
		      {
			mult_Minv(chi,ferm_conf,tp,iflav,mp->residue,eta);
			
			for(int iop=0;iop<nops;iop++)
			  summ_tens_dens(tens_dens[ind_op_flav(iop,iflav)],chiop,temp,ferm_conf,tp->backfield[iflav],shift[iop],mask[iop],eta);
		      }
		  }
		
		//plaquette and local charge
		double plaq=global_plaquette_lx_conf(smoothed_conf);
		local_topological_charge(topo_dens,smoothed_conf);
		
		for(int iop=0;iop<nops;iop++)
		  for(int iflav=0;iflav<nflavs;iflav++)
		    {
		      int iop_flav=ind_op_flav(iop,iflav);
		      
		      //final normalization and collapse
		      double_vector_prod_double((double*)(tens_dens[iop_flav]),(double*)(tens_dens[iop_flav]),1.0/nhits,loc_vol*2);
		      complex_vector_glb_collapse(tens[iop_flav],tens_dens[iop_flav],loc_vol);
		      //compute correlation with topocharge
		      compute_tens_dens_topo_correlation(spinpol_dens[iop_flav],tens_dens[iop_flav],topo_dens);
		      complex_vector_glb_collapse(spinpol[iop_flav],spinpol_dens[iop_flav],loc_vol);
		      
		      //total charge
		      double tot_charge;
		      double_vector_glb_collapse(&tot_charge,topo_dens,loc_vol);
		      double tot_charge2=double_vector_glb_norm2(topo_dens,1);
		      
		      fprintf(fout, "%d\t%d\t%d\t%d\t%d\t%d,%d\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\n",
			      iconf,icopy,imeas*meas_each,iconf,iflav,mp->operators[iop].first,mp->operators[iop].second,plaq,tot_charge,tot_charge2,
			      spinpol[iop_flav][RE],spinpol[iop_flav][IM],tens[iop_flav][RE],tens[iop_flav][IM]);
		    }
	      }
	  }
      }
  
    //free
    nissa_free(topo_dens);
    for(int iflav_op=0;iflav_op<nflavs*nops;iflav_op++)
      {
	nissa_free(spinpol_dens[iflav_op]);
	nissa_free(tens_dens[iflav_op]);
      }
    for(int eo=0;eo<2;eo++)
      {
	nissa_free(eta[eo]);
	nissa_free(chi[eo]);
	nissa_free(chiop[eo]);
	for(int itemp=0;itemp<2;itemp++)
	  nissa_free(temp[itemp][eo]);
	nissa_free(ferm_conf[eo]);
      }
    nissa_free(smoothed_conf);
    for(int is=0;is<nmeas*nhits;is++) nissa_free(source[is]);
    
    //close
    close_file(fout);
  }
  THREADABLE_FUNCTION_END
  
  std::string spinpol_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasSpinPol\n";
    os<<base_fermionic_meas_t::get_str(full);
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
