#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "geometry/geometry_mix.hpp"
#include "hmc/theory_pars.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "geometry/geometry_eo.hpp"

#include "spinpol.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  using namespace stag;	//In namespace stag there are all functions of mesons.cpp (L)
  typedef std::vector<std::pair<int,int> > op_list_t;
  
  THREADABLE_FUNCTION_7ARG(compute_tensorial_density, complex*,dens, complex**,loc_dens, theory_pars_t*,tp, quad_su3 **,conf, op_list_t,ops, int,nhits, double,residue)
  {
    GET_THREAD_ID();	//Define variable thread_id, say to the core of cpu how divide their work (L)
    int nop=ops.size();	//I changed nops in nop, because I've used nop later (L)
    
    //allocate noise and solution
    //Maybe this were respectivly source and sol (L). I didn't write them
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    
    //At the beginning define and alloc some quantity that I'll use later (L)
    int mask[nop],shift[nop];
    //Sol is the output of mult_Minv. temp is an internal variable of
    //apply_op, an input of apply_cov_shift. quark is the output of
    //apply_op (maybe I should chanche it). Source is eta (?). Then
    //allocate memory for all these variables. (L)
    color *source[2], *sol[2], *temp[2][2],*quark[nop][2];
    for(int eo=0;eo<2;eo++) source[eo]=nissa_malloc("source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) sol[eo]=nissa_malloc("sol",loc_volh+bord_volh,color);
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	temp[itemp][eo]=nissa_malloc("temp",loc_volh+bord_volh,color);
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	quark[iop][eo]=nissa_malloc("quark",loc_volh+bord_volh,color);
    
    //Create the mask before the hit loop. (L)
    for(int iop=0;iop<nop;iop++) {
      int spin=ops.at(iop).first;	//spin is index a of my note. Ops is an input vector of pair <a,b>. Same for taste (b). (L)
	int taste=ops.at(iop).second;
      shift[iop]=(spin^taste);
      mask[iop]=form_stag_op_pattern(spin,taste);	//This is the correct way to form the mask (?) (L)
      //if((shift[iop])&1) crash("operator %d (%d %d) has unmarched number of g0",iop,spin,taste);
      verbosity_lv3_master_printf(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
    }
    
    for(int iflav=0;iflav<tp->nflavs();iflav++)
      {
	if(tp->quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	
	//reset the local density
	for(int iop=0;iop<nop;iop++) vector_reset(loc_dens[iop+nop*iflav]);
	
	for(int ihit=0;ihit<nhits;ihit++) {
	    //We need this part to generate source eta (L)
	    //We don't need tso, see my notebook (L)
	    generate_fully_undiluted_eo_source(source,RND_Z4,-1);  //I called ori_source simply source, because I need only one (later) (L)
	   
	    //calculate sol (L)
	    mult_Minv(sol,conf,tp,iflav,residue,source);  //I modified residue (because of the input). (L)
	    
	    //Start with flav and op loop, calculate spinpol.(L)
	    for(int iop=0;iop<nop;iop++)
	    {
	      //calculate quark
	      apply_op(quark[iop],temp[0],temp[1],conf,tp->backfield[iflav],shift[iop],sol);
	      put_stag_phases(quark[iop],mask[iop]);  //From now use quark, but I have to link at loc_dens (L)
	    }
	    //Here, I have to put a routine that performs the sum on
	    //ihit. The result is put in loc_dens. DOn't know if this
	    //is the correct way to call complex_summassign
	    //(?). There's the problem of knowing how call quark, if
	    //with one or two indeces.(?) (L)
	    for(int iop=0;iop<nop;iop++)
	      for(int eo=0;eo<2;eo++)
		NISSA_PARALLEL_LOOP(ieo,0,loc_volh)	//loop on ieo and eo, i.e. on all eo-lattice (L)
		{
		  int ivol=loclx_of_loceo[eo][ieo];	//convert [eo][ieo] in lexicographically index ivol (L)
		  for(int ic=0;ic<NCOL;ic++)
		    //perform sum and scalar-product, see above (L)
		    //loc_dens every nop components changes flavour (L)
		    complex_summ_the_conj1_prod(loc_dens[iop+nop*iflav][ivol],source[eo][ieo][ic],quark[iop][eo][ieo][ic]);
		}
	    
	}
	//final normalization and collapse
	for(int iop=0;iop<nop;iop++)
	  {
	    double_vector_prod_double((double*)(loc_dens[iop+nop*iflav]),(double*)(loc_dens[iop+nop*iflav]),1.0/nhits,loc_vol*2);
	    complex_vector_glb_collapse(dens[iop+nop*iflav],loc_dens[iop+nop*iflav],loc_vol);
	  }
      }
    
    //free
    for(int par=0;par<2;par++)
      {
	nissa_free(rnd[par]);
	nissa_free(chi[par]);
      }
    //Deallocate memory used for my variables (L)
    for(int eo=0;eo<2;eo++) nissa_free(source[eo]);
    for(int eo=0;eo<2;eo++) nissa_free(sol[eo]);
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	nissa_free(temp[itemp][eo]);
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	nissa_free(quark[iop][eo]);
  }
  THREADABLE_FUNCTION_END
  
  //compute the spin-polarization for all flavors
  void measure_spinpol(quad_su3 **ignore_ferm_conf,theory_pars_t &tp,spinpol_meas_pars_t &mp,int iconf,int conf_created,stout_pars_t &stout_pars,quad_su3 **glu_conf)
  {
    verbosity_lv1_master_printf("Evaluating spinpol\n");
    
    //    if(mp.use_ferm_conf_for_gluons or glu_conf==NULL) glu_conf=ferm_conf;	//We don't use ferm_conf anymore (L)
    
    smooth_pars_t &sp=mp.smooth_pars;
    int nflavs=tp.nflavs();
    int nop=mp.nops();
    
    //open the file
    FILE *fout=open_file(mp.path,conf_created?"w":"a");
    
    //allocate point and local results
    double *topo_dens=nissa_malloc("topo_dens",loc_vol,double);
    complex *spinpol_dens=nissa_malloc("spinpol_dens",loc_vol,complex);
    complex tens[nflavs*nop];
    complex *tens_dens[nflavs*nop];	//This will be local_dens, when we call compute_tens_dens (L)
    for(int iflav_op=0;iflav_op<nflavs*nop;iflav_op++) tens_dens[iflav_op]=nissa_malloc("tens_dens",loc_vol+bord_vol,complex);
    
    //allocate the smoothed confs
    int nmeas=mp.smooth_pars.nmeas_nonzero()+1;
    quad_su3 *smoothed_conf[nmeas];
    quad_su3 *ferm_conf[nmeas][2];
    quad_su3 *gauge_conf=nissa_malloc("gauge_conf",loc_vol+bord_vol,quad_su3);
    for(int imeas=0;imeas<nmeas;imeas++)
      {
	smoothed_conf[imeas]=nissa_malloc(combine("smoothed_conf_%d",imeas).c_str(),loc_vol+bord_vol,quad_su3);
	for(int eo=0;eo<2;eo++) ferm_conf[imeas][eo]=nissa_malloc("ferm_conf",loc_volh+bord_volh+edge_volh,quad_su3);
      }
    if(nmeas==0) crash("nmeas cannot be 0");
    
    //smooth
    int imeas=0;
    int nsmooth=0;
    std::vector<int> nsmooth_meas(nmeas);
    bool finished;
    do
      {
	verbosity_lv2_master_printf("Meas: %d/%d, nsmooth: %d\n",imeas,nmeas,nsmooth);
	
	if(imeas==0)
	  {
	    paste_eo_parts_into_lx_vector(smoothed_conf[0],glu_conf);
	    finished=(nmeas==1);
	  }
	else
	  {
	    vector_copy(smoothed_conf[imeas],smoothed_conf[imeas-1]);
	    finished=smooth_lx_conf_until_next_meas(smoothed_conf[imeas],sp,nsmooth);
	  }
	
	split_lx_vector_into_eo_parts(ferm_conf[imeas],smoothed_conf[imeas]);
	stout_smear(ferm_conf[imeas],ferm_conf[imeas],&stout_pars);
	
	nsmooth_meas[imeas]=nsmooth;
	imeas++;
      }
    while(not finished);
    
    //compute the topological charge and the product of topological and tensorial density
    int ncopies=mp.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      for(int imeas=0;imeas<nmeas;imeas++)
	{
	  verbosity_lv1_master_printf("Computing copy %d/%d smooth %d/%d, t %d/%d\n",icopy,ncopies,imeas,nmeas,nsmooth_meas[imeas],nsmooth);
	  
	  //evaluate the tensorial density for all quarks
	  compute_tensorial_density(tens,tens_dens,&tp,ferm_conf[imeas],mp.operators,mp.nhits,mp.residue);
	  
	  for(int igauge_conf=0;igauge_conf<2;igauge_conf++)
	    {
	      //take the gauge conf, either the smoothed conf or the fermionic conf
	      if(igauge_conf==0) vector_copy(gauge_conf,smoothed_conf[imeas]);
	      else               paste_eo_parts_into_lx_vector(gauge_conf,ferm_conf[imeas]);
	      
	      //plaquette and local charge
	      double plaq=global_plaquette_lx_conf(gauge_conf);
	      local_topological_charge(topo_dens,gauge_conf);
	      //total charge
	      double tot_charge;
	      double_vector_glb_collapse(&tot_charge,topo_dens,loc_vol);
	      double tot_charge2=double_vector_glb_norm2(topo_dens,1);
	      
	      //topo-tens
	      for(int iop=0;iop<mp.nops();iop++)
		for(int iflav=0;iflav<tp.nflavs();iflav++)
		  {
		    GET_THREAD_ID();
		    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) complex_prod_double(spinpol_dens[ivol],tens_dens[iop+nop*iflav][ivol],topo_dens[ivol]);
		    THREAD_BARRIER();
		    
		    complex spinpol;
		    complex_vector_glb_collapse(spinpol,spinpol_dens,loc_vol);
		    master_fprintf(fout, "%d\t%d\t%d\t%d\t%d\t%d,%d\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\n",
				   iconf,icopy,nsmooth_meas[imeas],igauge_conf,iflav,mp.operators[iop].first,mp.operators[iop].second,plaq,tot_charge,tot_charge2,
				   spinpol[RE],spinpol[IM],tens[iop+nop*iflav][RE],tens[iop+nop*iflav][IM]);
		  }
	    }
	}
    
    //free
    nissa_free(topo_dens);
    for(int iflav_op=0;iflav_op<nflavs*nop;iflav_op++) nissa_free(tens_dens[iflav_op]);
    nissa_free(spinpol_dens);
    for(int imeas=0;imeas<nmeas;imeas++)
      {
	nissa_free(smoothed_conf[imeas]);
	for(int eo=0;eo<2;eo++) nissa_free(ferm_conf[imeas][eo]);
      }
    nissa_free(gauge_conf);
    
    //close
    close_file(fout);
  }
  
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
