#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "linalgs/reduce.hpp"
#include "mesons.hpp"
#include "new_types/su3.hpp"
#include "operations/gauge_fixing.hpp"
#include "routines/mpi_routines.hpp"
#include "stag.hpp"

namespace nissa {
using namespace stag;

namespace
{
  int nop;
  int ncombo;
  int nflavs;
}
  
  //compute the index where to store
  inline int icombo(int iflav,int iop_so,int iop_si,const GlbCoord& t)
  {
    return (t+glbSize(tDir)*(iop_si+nop*(iop_so+nop*iflav)))();
  }
  
  //compute correlation functions for staggered mesons, arbitary taste and spin
  void compute_meson_corr(complex* corr,eo_ptr<quad_su3> conf,theory_pars_t* tp,meson_corr_meas_pars_t* meas_pars)
  {
eo_ptr<color> ori_source,source,sol,quark[nop],temp[2];
    for(int eo=0;eo<2;eo++) ori_source[eo]=nissa_malloc("ori_source",locVolhWithBord.nastyConvert(),color);
    for(int eo=0;eo<2;eo++) source[eo]=nissa_malloc("source",locVolhWithBord.nastyConvert(),color);
    for(int eo=0;eo<2;eo++) sol[eo]=nissa_malloc("sol",locVolhWithBord.nastyConvert(),color);
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	quark[iop][eo]=nissa_malloc("quark",locVolhWithBord.nastyConvert(),color);
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	temp[itemp][eo]=nissa_malloc("temp",locVolhWithBord.nastyConvert(),color);
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
	//if((shift[iop])&1) crash("operator %d (%d %d) has unmatched number of g0",iop,spin,taste);
	verbosity_lv3_master_printf(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
      }
    
    vector_reset(corr);
    
    //measure the putpourri for each quark
    for(int ihit=0;ihit<meas_pars->nhits;ihit++)
      {
	verbosity_lv2_master_printf("Computing hit %d/%d\n",ihit,meas_pars->nhits);
	
	//generate tso
	int source_coord;
	source_coord=rnd_get_unif(&glb_rnd_gen,0,glbTimeSize());
	verbosity_lv2_master_printf("tsource: %d\n",source_coord);
	
	//generate source
	generate_fully_undiluted_eo_source(ori_source,meas_pars->rnd_type,source_coord);
	
	// for(int iflav=0;iflav<nflavs;iflav++)
	//   {
	//     for(int iop=0;iop<nop;iop++)
	//       {
	// 	apply_shift_op(source,temp[0],temp[1],conf,tp->backfield[iflav],shift[iop],ori_source);
	// 	put_stag_phases(source,mask[iop]);
	// 	mult_Minv(quark[iop],conf,tp,iflav,meas_pars->residue,source);
	//       }
	    
	//     /// Sink
	//     for(int iop=0;iop<nop;iop++)
	//       {
	// 	apply_shift_op(quark0s[iop],temp[0],temp[1],conf,tp->backfield[iflav],shift[iop],quark[0]);
	// 	put_stag_phases(quark0s[iop],mask[iop]);
	//       }
	    
	//     for(int iop_so=0;iop_so<nop;iop_so++)
	//       for(int iop_si=0;iop_si<nop;iop_si++)
	// 	{
	// 	  complex *loc_contr=get_reducing_buffer<complex>(locVol);
		  
	// 	  for(int eo=0;eo<2;eo++)
	// 	    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	// 	      {
	// 		// int ivol=loclx_of_loceo[eo][ieo];
	// 		// int t=(glb_coord_of_loclx[ivol][0]-tso+glb_size[0])%glb_size[0];
	// 		// for(int ic=0;ic<NCOL;ic++)
	// 		//   complex_summ_the_conj1_prod(loc_corr[icombo(iflav,iop,t)],quark[0][eo][ieo][ic],quark[iop][eo][ieo][ic]);
	// 	      }
	// 	NISSA_PARALLEL_LOOP_END;
	// 	THREAD_BARRIER();
	//       }
	//   }
	
	//reduce
	crash("#warning reimplement glb_threads_reduce_double_vect((double*)loc_corr,2*ncombo");
	crash("#warning reimplement if(IS_MASTER_THREAD) glb_nodes_reduce_complex_vect(corr,loc_corr,ncombo");
      }
    
    // for(int eo=0;eo<2;eo++)
    //   {
    // 	nissa_free(ori_source[eo]);
    // 	nissa_free(source[eo]);
    // 	nissa_free(sol[eo]);
    //   }
    // for(int iop=0;iop<nop;iop++)
    //   {
    // 	for(int eo=0;eo<2;eo++)
    // 	  nissa_free(quark[iop][eo]);
    // 	for(int eo=0;eo<2;eo++)
    // 	  nissa_free(quark0s[iop][eo]);
    //   }
    // for(int itemp=0;itemp<2;itemp++)
    //   for(int eo=0;eo<2;eo++)
    // 	nissa_free(temp[itemp][eo]);
    
    // delete[] quark;
    // delete[] quark0s;
  }
  
  //compute and print
  void measure_meson_corr(eo_ptr<quad_su3> ext_conf,theory_pars_t &tp,meson_corr_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    nop=meas_pars.mesons.size();
    nflavs=tp.nflavs();
    ncombo=icombo(nflavs-1,nop-1,nop-1,glbTimeSize-1)+1;
    double norm=1.0/(meas_pars.nhits*glbSpatVol());
    complex *corr=nissa_malloc("corr",ncombo,complex);
    
    //measure the meson corrs for each quark
    int ncopies=meas_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
	verbosity_lv2_master_printf("Computing copy %d/%d\n",icopy,ncopies);
	
	compute_meson_corr(corr,ext_conf,&tp,&meas_pars);
	
	//open the file, allocate point result and source
	FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
	for(int iop_so=0;iop_so<nop;iop_so++)
	  for(int iop_si=0;iop_si<nop;iop_si++)
	    for(int iflav=0;iflav<nflavs;iflav++)
	      {
		int spin_so=meas_pars.mesons[iop_so].first;
		int taste_so=meas_pars.mesons[iop_so].second;
		int spin_si=meas_pars.mesons[iop_si].first;
		int taste_si=meas_pars.mesons[iop_si].second;
		master_fprintf(file," # conf %d ;"
			       " iop_so %d , spin_so %d , taste_so %d ;"
			       " iop_si %d , spin_si %d , taste_si %d ;"
			       " flv = %d , m = %lg\n",
			       iconf,iop_so,spin_so,taste_so,iop_si,spin_si,taste_si,iflav,tp.quarks[iflav].mass);
		for(GlbCoord t=0;t<glbTimeSize;t++)
		  {
		    int ic=icombo(iflav,iop_so,iop_si,t);
		    master_fprintf(file,"%d %+16.16lg %+16.16lg\n",t,corr[ic][RE]*norm,corr[ic][IM]*norm);
		  }
		master_fprintf(file,"\n");
	      }
	close_file(file);
	
	nissa_free(corr);
      }
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
