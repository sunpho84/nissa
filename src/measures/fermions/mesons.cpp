#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "mesons.hpp"
#include "new_types/su3.hpp"
#include "operations/gauge_fixing.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/thread.hpp"
#include "stag.hpp"

namespace nissa {
using namespace stag;

namespace {
int nop;
int ncombo;
int nflavs;
  }
  
  //compute the index where to store
  inline int icombo(int iflav,int iop,int t)
  {return t+glb_size[0]*(iop+nop*iflav);}
  
  //compute correlation functions for staggered mesons, arbitary taste and spin
  THREADABLE_FUNCTION_4ARG(compute_meson_corr, complex*,corr, quad_su3**,conf, theory_pars_t*,tp, meson_corr_meas_pars_t*,meas_pars)
  {
    GET_THREAD_ID();
    
    //allocate
    color *ori_source[2],*source[2],*sol[2],*quark[nop][2],*temp[2][2];
    for(int eo=0;eo<2;eo++) ori_source[eo]=nissa_malloc("ori_source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) source[eo]=nissa_malloc("source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) sol[eo]=nissa_malloc("sol",loc_volh+bord_volh,color);
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	quark[iop][eo]=nissa_malloc("quark",loc_volh+bord_volh,color);
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	temp[itemp][eo]=nissa_malloc("temp",loc_volh+bord_volh,color);
    complex *loc_corr=new complex[ncombo];
    memset(loc_corr,0,sizeof(complex)*ncombo);
    
    //form the masks
    int mask[nop],shift[nop];
    for(int iop=0;iop<nop;iop++)
      {
	int spin=meas_pars->mesons[iop].first;
	int taste=meas_pars->mesons[iop].second;
	shift[iop]=(spin^taste);
	mask[iop]=form_stag_meson_pattern_with_g5g5(spin,taste);
	//if((shift[iop])&1) crash("operator %d (%d %d) has unmarched number of g0",iop,spin,taste);
	verbosity_lv3_master_printf(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
      }
    
    //measure the putpourri for each quark
    int ncopies=meas_pars->ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
	for(int ihit=0;ihit<meas_pars->nhits;ihit++)
	  {
	    verbosity_lv2_master_printf("Computing copy %d/%d hit %d/%d\n",icopy,ncopies,ihit,meas_pars->nhits);
	    
	    //generate tso
	    int tso;
	    if(IS_MASTER_THREAD) tso=rnd_get_unif(&glb_rnd_gen,0,glb_size[0]);
	    THREAD_BROADCAST(tso,tso);
	    verbosity_lv2_master_printf("tsource: %d\n",tso);
	    
	    //generate source
	    generate_fully_undiluted_eo_source(ori_source,meas_pars->rnd_type,tso);
	    
	    for(int iflav=0;iflav<nflavs;iflav++)
	      {
		for(int iop=0;iop<nop;iop++)
		  {
		    apply_shift_op(source,temp[0],temp[1],conf,tp->backfield[iflav],shift[iop],ori_source);
		    put_stag_phases(source,mask[iop]);
		    mult_Minv(sol,conf,tp,iflav,meas_pars->residue,source);
		    apply_shift_op(quark[iop],temp[0],temp[1],conf,tp->backfield[iflav],shift[iop],sol);
		    put_stag_phases(quark[iop],mask[iop]);
		  }
		
		for(int iop=0;iop<nop;iop++)
		  for(int eo=0;eo<2;eo++)
		    NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
		      {
			int ivol=loclx_of_loceo[eo][ieo];
			int t=(glb_coord_of_loclx[ivol][0]-tso+glb_size[0])%glb_size[0];
			for(int ic=0;ic<NCOL;ic++)
			  complex_summ_the_conj1_prod(loc_corr[icombo(iflav,iop,t)],quark[0][eo][ieo][ic],quark[iop][eo][ieo][ic]);
		      }
		THREAD_BARRIER();
	      }
	  }
	
	//reduce
	glb_threads_reduce_double_vect((double*)loc_corr,2*ncombo);
	if(IS_MASTER_THREAD) glb_nodes_reduce_complex_vect(corr,loc_corr,ncombo);
      }
    
    for(int eo=0;eo<2;eo++)
      {
	nissa_free(ori_source[eo]);
	nissa_free(source[eo]);
	nissa_free(sol[eo]);
      }
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	nissa_free(quark[iop][eo]);
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	nissa_free(temp[itemp][eo]);
    delete[] loc_corr;
  }
  THREADABLE_FUNCTION_END
  
  //compute and print
  void measure_meson_corr(quad_su3 **ext_conf,theory_pars_t &tp,meson_corr_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    nop=meas_pars.mesons.size();
    nflavs=tp.nflavs();
    ncombo=icombo(nflavs-1,nop-1,glb_size[0]-1)+1;
    double norm=1.0/(meas_pars.nhits*glb_spat_vol);
    complex *corr=nissa_malloc("corr",ncombo,complex);
    
    compute_meson_corr(corr,ext_conf,&tp,&meas_pars);
    
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    for(int iop=0;iop<nop;iop++)
      for(int iflav=0;iflav<nflavs;iflav++)
	{
	  int spin=meas_pars.mesons[iop].first;
	  int taste=meas_pars.mesons[iop].second;
	  master_fprintf(file," # conf %d ; iop %d , spin %d , taste %d ; flv = %d , m = %lg\n",
			 iconf,iop,spin,taste,iflav,tp.quarks[iflav].mass);
	  for(int t=0;t<glb_size[0];t++)
	    {
	      int ic=icombo(iflav,iop,t);
	      master_fprintf(file,"%d %+16.16lg %+16.16lg\n",t,corr[ic][RE]*norm,corr[ic][IM]*norm);
	    }
	    master_fprintf(file,"\n");
	  }
    close_file(file);
    
    nissa_free(corr);
  }
  
  //nucleon correlators
  std::string meson_corr_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMesonCorrs\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(mesons.size() or full)
      {
	os<<" Operators\t=\t{";
	for(size_t i=0;i<mesons.size();i++)
	  {
	    os<<"("<<mesons[i].first<<","<<mesons[i].second<<")";
	    if(i!=mesons.size()-1) os<<",";
	  }
	os<<"}\n";
      }
    
    return os.str();
  }
  } // namespace nissa
