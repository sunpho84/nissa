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

#include "spinpol.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  typedef std::vector<std::pair<int,int> > op_list_t;
  
  THREADABLE_FUNCTION_7ARG(compute_tensorial_density, complex*,dens, complex**,loc_dens, theory_pars_t*,tp, quad_su3 **,conf, op_list_t,ops, int,nhits, double,residue)
  {
    int nops=ops.size();
    
    //allocate noise and solution
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    
    for(int iflav=0;iflav<tp->nflavs();iflav++)
      {
	if(tp->quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	
	//reset the local density
	for(int iop=0;iop<nops;iop++) vector_reset(loc_dens[iop+nops*iflav]);
	for(int ihit=0;ihit<nhits;ihit++)
	  {
	    //AOH QUA C'E' DA SCRIVERE
	    
	  }
	//final normalization and collapse
	for(int iop=0;iop<nops;iop++)
	  {
	    double_vector_prod_double((double*)(loc_dens[iop+nops*iflav]),(double*)(loc_dens[iop+nops*iflav]),1.0/nhits,loc_vol*2);
	    complex_vector_glb_collapse(dens[iop+nops*iflav],loc_dens[iop+nops*iflav],loc_vol);
	  }
      }
    
    //free
    for(int par=0;par<2;par++)
      {
	nissa_free(rnd[par]);
	nissa_free(chi[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //compute the spin-polarization for all flavors
  void measure_spinpol(quad_su3 **ferm_conf,theory_pars_t &tp,spinpol_meas_pars_t &mp,int iconf,int conf_created,quad_su3 **glu_conf)
  {
    if(mp.use_ferm_conf_for_gluons or glu_conf==NULL) glu_conf=ferm_conf;
    
    smooth_pars_t &sp=mp.smooth_pars;
    int nflavs=tp.nflavs();
    int nops=mp.nops();
    
    //open the file
    FILE *fout=open_file(mp.path,conf_created?"w":"a");
    
    //allocate point and local results
    double *topo_dens=nissa_malloc("topo_dens",loc_vol,double);
    complex *spinpol_dens=nissa_malloc("spinpol_dens",loc_vol,complex);
    complex tens[nflavs*nops];
    complex *tens_dens[nflavs*nops];
    for(int iflav_op=0;iflav_op<nflavs*nops;iflav_op++) tens_dens[iflav_op]=nissa_malloc("tens_dens",loc_vol+bord_vol,complex);
    
    //evaluate the tensorial density for all quarks
    compute_tensorial_density(tens,tens_dens,&tp,ferm_conf,mp.operators,mp.nhits,mp.residue);
    
    //compute the topological charge and the product of topological and tensorial density
    quad_su3 *smoothed_conf=nissa_malloc("smoothed_conf",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_vector(smoothed_conf,glu_conf);
    double t=0,tnext_meas=sp.meas_each;
    bool finished;
    do
      {
	//plaquette and local charge
	double plaq=global_plaquette_lx_conf(smoothed_conf);
	local_topological_charge(topo_dens,smoothed_conf);
	//total charge
	double tot_charge;
	double_vector_glb_collapse(&tot_charge,topo_dens,loc_vol);
	
	//topo-tens
	for(int iop=0;iop<mp.nops();iop++)
	  for(int iflav=0;iflav<tp.nflavs();iflav++)
	    {
	      GET_THREAD_ID();
	      NISSA_PARALLEL_LOOP(ivol,0,loc_vol) complex_prod_double(spinpol_dens[ivol],tens_dens[iop+nops*iflav][ivol],topo_dens[ivol]);
	      THREAD_BARRIER();
	      
	      complex spinpol;
	      complex_vector_glb_collapse(spinpol,spinpol_dens,loc_vol);
	      master_fprintf(fout, "%d\t%lg\t%d\t%d,%d\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\n",iconf,t,iflav,mp.operators[iop].first,mp.operators[iop].second,plaq,tot_charge,spinpol[RE],spinpol[IM]);
	    }
	finished=smooth_lx_conf_until_next_meas(smoothed_conf,sp,t,tnext_meas);
      }
    while(!finished);
    
    //free
    nissa_free(topo_dens);
    for(int iflav_op=0;iflav_op<nflavs*nops;iflav_op++) nissa_free(tens_dens[iflav_op]);
    nissa_free(spinpol_dens);
    nissa_free(smoothed_conf);
    
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
