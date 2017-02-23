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
  THREADABLE_FUNCTION_7ARG(compute_tensorial_density, complex*,dens, complex**,loc_dens, theory_pars_t*,tp, quad_su3 **,conf, std::vector<int>,dirs, int,nhits, double,residue)
  {
    //allocate noise and solution
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    
    for(int iflav=0;iflav<tp->nflavs();iflav++)
      {
	if(tp->quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	
	//reset the local density
	vector_reset(loc_dens[iflav]);
	for(int ihit=0;ihit<nhits;ihit++)
	  {
	  }
	//final normalization and collapse
	double_vector_prod_double((double*)(loc_dens[iflav]),(double*)(loc_dens[iflav]),1.0/nhits,loc_vol*2);
	complex_vector_glb_collapse(dens[iflav],loc_dens[iflav],loc_vol);
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
  void measure_spinpol(quad_su3 **ferm_conf,quad_su3 **glu_conf,theory_pars_t &tp,spinpol_meas_pars_t &mp,int iconf,int conf_created)
  {
    if(mp.use_ferm_conf_for_gluons) glu_conf=ferm_conf;
    
    smooth_pars_t &sp=mp.smooth_pars;
    int nflavs=tp.nflavs();
    int ndirs=mp.ndirs();
    
    //open the file
    FILE *fout=open_file(mp.path,conf_created?"w":"a");
    
    //allocate point and local results
    double *topo_dens=nissa_malloc("topo_dens",loc_vol,double);
    complex *spinpol_dens=nissa_malloc("spinpol_dens",loc_vol,complex);
    complex tens[nflavs*ndirs];
    complex *tens_dens[nflavs*ndirs];
    for(int iflav_dir=0;iflav_dir<nflavs*ndirs;iflav_dir++) tens_dens[iflav_dir]=nissa_malloc("tens_dens",loc_vol+bord_vol,complex);
    
    //evaluate the tensorial density for all quarks
    compute_tensorial_density(tens,tens_dens,&tp,ferm_conf,mp.dirs,mp.nhits,mp.residue);
    
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
	for(int idir=0;idir<mp.ndirs();idir++)
	  for(int iflav=0;iflav<tp.nflavs();iflav++)
	    {
	      GET_THREAD_ID();
	      NISSA_PARALLEL_LOOP(ivol,0,loc_vol) complex_prod_double(spinpol_dens[ivol],tens_dens[idir+ndirs*iflav][ivol],topo_dens[ivol]);
	      THREAD_BARRIER();
	      
	      complex spinpol;
	      complex_vector_glb_collapse(spinpol,spinpol_dens,loc_vol);
	      master_fprintf(fout, "%d\t%lg\t%d\t%d\t%+16.16lg\t%+16.16lg\t%+16.16lg\t%+16.16lg\n",iconf,t,iflav,idir,plaq,tot_charge,spinpol[RE],spinpol[IM]);
	    }
	finished=smooth_lx_conf_until_next_meas(smoothed_conf,sp,t,tnext_meas);
      }
    while(!finished);
    
    //free
    nissa_free(topo_dens);
    for(int iflav_dir=0;iflav_dir<nflavs*ndirs;iflav_dir++) nissa_free(tens_dens[iflav_dir]);
    nissa_free(spinpol_dens);
    nissa_free(smoothed_conf);
    
    //close
    close_file(fout);
  }
  
  std::string spinpol_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasTop\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(dirs.size())
      {
	os<<"Dirs\t=\t{"<<dirs[0];
	for(size_t idir=1;idir<dirs.size();idir++) os<<","<<dirs[idir];
	os<<"}";
      }
    os<<smooth_pars.get_str(full);
    
    return os.str();
  }
}
