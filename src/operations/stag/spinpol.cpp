#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  THREADABLE_FUNCTION_7ARG(compute_tensorial_density, complex*,dens, complex**,loc_dens, theory_pars_t*,tp, quad_su3 **,conf, int,dir, int,nhits, double,residue)
  {
    //allocate noise and solution
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    
    for(int iflav=0;iflav<tp->nflavs;iflav++)
      {
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
  void measure_spinpol(quad_su3 **ferm_conf,quad_su3 **glu_conf,theory_pars_t &tp,int iconf,int conf_created)
  {
    spinpol_meas_pars_t *sp=&tp.spinpol_meas_pars;
    if(sp->use_ferm_conf_for_gluons) glu_conf=ferm_conf;
    
    //count the number of cooled levels
    cool_pars_t *cp=&sp->cool_pars;
    int ncool_meas=cp->nsteps/cp->meas_each+1;
    verbosity_lv3_master_printf("ncool_meas: %d\n",ncool_meas);
    
    //allocate point and local results
    double *topo_dens=nissa_malloc("topo_dens",loc_vol,double);
    complex tens[tp.nflavs];
    complex *tens_dens[tp.nflavs];
    complex *spinpol_dens=nissa_malloc("spinpol_dens",loc_vol,complex);
    for(int iflav=0;iflav<tp.nflavs;iflav++) tens_dens[iflav]=nissa_malloc("tens_dens_loc",loc_vol+bord_vol,complex);
    
    //evaluate the tensorial density for all quarks
    compute_tensorial_density(tens,tens_dens,&tp,ferm_conf,sp->dir,sp->nhits,sp->residue);
    
    //compute the topological charge and the product of topological and tensorial density
    quad_su3 *cooled_conf=nissa_malloc("cooled_conf",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_conf(cooled_conf,glu_conf);
    double topo[ncool_meas];
    complex spinpol[ncool_meas][tp.nflavs];
    for(int icool=0;icool<=cp->nsteps;icool++)
      {
	if(icool%cp->meas_each)
	  {
	    int imeas=icool/cp->meas_each;
	    
	    //topological charge
	    local_topological_charge(topo_dens,cooled_conf);
	    double_vector_glb_collapse(topo+imeas,topo_dens,loc_vol);
	    
	    //topo-tens
	    for(int iflav=0;iflav<tp.nflavs;iflav++)
	      {
		GET_THREAD_ID();
		NISSA_PARALLEL_LOOP(ivol,0,loc_vol) complex_prod_double(spinpol_dens[ivol],tens_dens[iflav][ivol],topo_dens[ivol]);
		THREAD_BARRIER();
		
		complex_vector_glb_collapse(spinpol[imeas][iflav],spinpol_dens,loc_vol);
	      }
	  }
	
	//cool if needed
	if(icool!=cp->nsteps) cool_lx_conf(cooled_conf,cp->gauge_action,cp->overrelax_flag,cp->overrelax_exp);
      }
    
    //free
    for(int icool=0;icool<ncool_meas;icool++) nissa_free(topo_dens);
    for(int iflav=0;iflav<tp.nflavs;iflav++) nissa_free(tens_dens[iflav]);
    nissa_free(spinpol_dens);
    nissa_free(cooled_conf);
    
    //////////////////////////////////// output //////////////////////////////////
    
    //open the file and write header
    FILE *fout=open_file(tp.spinpol_meas_pars.path,conf_created?"w":"a");
    master_fprintf(fout," # conf %d, nhits %d, dir %d, residue %lg\n\n",iconf,sp->nhits,sp->dir,sp->residue);
    
    //write tensorial density alone
    for(int iflav=0;iflav<tp.nflavs;iflav++) master_fprintf(fout," tdens flav %d:\t%+016.016lg\t%+016.016lg\n",iflav,tens[iflav][RE],tens[iflav][IM]);
    master_fprintf(fout,"\n");
    
    //write topological charge density
    for(int imeas=0;imeas<ncool_meas;imeas++) master_fprintf(fout," topocharge ncool %d:\t%+016.016lg\n",imeas*cp->meas_each,topo[imeas]);
    master_fprintf(fout,"\n");
    
    //write tensorial*topo density
    for(int imeas=0;imeas<ncool_meas;imeas++)
      for(int iflav=0;iflav<tp.nflavs;iflav++)
	master_fprintf(fout," spinpol flav %d cool %d:\t%+016.016lg\t%+016.016lg\n",iflav,imeas*cp->meas_each,spinpol[imeas][iflav][RE],spinpol[imeas][iflav][IM]);
    master_fprintf(fout,"\n\n\n");
    
    //close
    close_file(fout);
  }
}
