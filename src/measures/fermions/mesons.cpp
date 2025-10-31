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
#include "linalgs/reduce.hpp"
#include "mesons.hpp"
#include "new_types/su3.hpp"
#include "operations/gauge_fixing.hpp"
#include "routines/mpi_routines.hpp"
#include "stag.hpp"

namespace nissa
{
  using namespace stag;
  
  namespace
  {
    int nop;
    int ncombo;
    int nflavs;
  }
  
  //compute the index where to store
  inline int icombo(const int& iflav_so,
		    const int& iflav_si,
		    const int& iop_so,
		    const int& iop_si,
		    const int& i_dir,
		    const int& dir)
  {
    return i_dir
		+glbSize[dir]*(iop_si
			       +nop*(iop_so
				     +nop*(iflav_si
					   +nflavs*iflav_so)));
  }
  
  void compute_meson_corr(complex* corr,
			  const EoField<quad_su3>& conf,
			  const theory_pars_t& tp,
			  const meson_corr_meas_pars_t& meas_pars)
  {
    const int& dir=
      meas_pars.dir;
    
    //allocate
    EoField<color> ori_source("ori_source",WITH_HALO);
    EoField<color> source("source",WITH_HALO);
    EoField<color> sol("sol",WITH_HALO);
    std::vector<EoField<color>> temp(2,{"temp",WITH_HALO});
    std::vector<EoField<color>> quark(nflavs*nop,{"quark",WITH_HALO});
    std::vector<EoField<color>> quark0s(nflavs*nop,{"quark0s",WITH_HALO});
    
    //form the masks
    int mask[nop],shift[nop];
    for(int iop=0;iop<nop;iop++)
      {
	const auto& [spin,taste]=
	  meas_pars.mesons[iop];
	
	shift[iop]=(spin^taste);
	mask[iop]=form_stag_meson_pattern_with_g5g5(spin,taste);
	//if((shift[iop])&1) CRASH("operator %d (%d %d) has unmatched number of g0",iop,spin,taste);
	VERBOSITY_LV3_MASTER_PRINTF(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
      }
    
    for(int ic=0;ic<ncombo;ic++)
      complex_put_to_zero(corr[ic]);
    
    //measure the putpourri for each quark
    for(int ihit=0;ihit<meas_pars.nhits;ihit++)
      {
	VERBOSITY_LV2_MASTER_PRINTF("Computing hit %d/%d\n",ihit,meas_pars.nhits);
	
	//generate tso
	const int source_coord=
	  rnd_get_unif(&glb_rnd_gen,0,glbSize[dir]);
	VERBOSITY_LV2_MASTER_PRINTF("source: %d\n",source_coord);
	
	//generate source
	generate_fully_undiluted_eo_source(ori_source,meas_pars.rnd_type,source_coord,dir);
	
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    for(int iop=0;iop<nop;iop++)
	      {
			const int idx=iflav*nop+iop;
			apply_shift_op(source,temp[0],temp[1],conf,tp.backfield[iflav],shift[iop],ori_source);
			put_stag_phases(source,mask[iop]);
			mult_Minv(quark[idx],conf,tp,iflav,meas_pars.residue,source);
	      }
	  }
	    /// Sink
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    for(int iop=0;iop<nop;iop++)
	      {
			const int idx = iflav*nop + iop;
			apply_shift_op(quark0s[idx],temp[0],temp[1],conf,tp.backfield[iflav],shift[iop],quark[idx]);
			put_stag_phases(quark0s[idx],mask[iop]);
	      }
	  }
	
	// contract all flavs and ops
	for(int iflav_so=0;iflav_so<nflavs;iflav_so++)
	  for(int iflav_si=0;iflav_si<nflavs;iflav_si++)
	    for(int iop_so=0;iop_so<nop;iop_so++)
	      for(int iop_si=0;iop_si<nop;iop_si++)
		{
		  const int idx_so=iflav_so*nop+iop_so;
		  const int idx_si=iflav_si*nop+iop_si;
		  
		  LxField<complex> loc_contr("loc_contr");
		  
		  for(int eo=0;eo<2;eo++)
		    {
		      PAR(0,
			  locVolh,
			  CAPTURE(eo,
				  q0s=quark0s[idx_si].getReadable(),
				  q=quark[idx_so].getReadable(),
				  TO_WRITE(loc_contr)),
			  ieo,
			  {
			    const int ivol=
			      loclx_of_loceo[eo][ieo];
			    
			    color_scalar_prod(loc_contr[ivol],
					      q0s[eo][ieo],
					      q[eo][ieo]);
			  });
		    }
		  
		  complex unshiftedGlbContr[glbSize[dir]];
		  glb_reduce(unshiftedGlbContr,loc_contr,locVol,glbSize[dir],locSize[dir],glbCoordOfLoclx[0][dir]);
		  
		  for(int glb_idir=0;glb_idir<glbSize[dir];glb_idir++)
		    {
		      /// Distance from source
		      const int d_dir=
			(glb_idir-source_coord+glbSize[dir])%glbSize[dir];
		      
		      complex_summassign(corr[icombo(iflav_so,iflav_si,iop_so,iop_si,d_dir,dir)],unshiftedGlbContr[glb_idir]);
		    }
		}
      }
    
  }
  
  /// Compute and print
  void measure_meson_corr(const EoField<quad_su3>& ext_conf,
			  const theory_pars_t& tp,
			  const meson_corr_meas_pars_t& meas_pars,
			  const int& iconf,
			  const int& conf_created)
  {
    nop=meas_pars.mesons.size();
    nflavs=tp.nflavs();
    
    const int& dir=meas_pars.dir;
    
    VERBOSITY_LV1_MASTER_PRINTF("Meson correlator: type=%s, dir=%d, extent=%d\n",(dir==0)?"temporal":"spatial",dir,glbSize[dir]);
    
    ncombo=icombo(nflavs-1,nflavs-1,nop-1,nop-1,glbSize[dir]-1,dir)+1;
    const double orthVol=double(glbVol)/double(glbSize[dir]);
    double norm=1.0/(meas_pars.nhits*orthVol);
    
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    complex corr[ncombo];
    
    //measure the meson corrs for each quark
    for(int ncopies=meas_pars.ncopies,icopy=0;icopy<ncopies;icopy++)
      {
	VERBOSITY_LV2_MASTER_PRINTF("Computing copy %d/%d\n",icopy,ncopies);
	
	compute_meson_corr(corr,ext_conf,tp,meas_pars);
	
	//open the file, allocate point result and source
	for(int iop_so=0;iop_so<nop;iop_so++)
	  for(int iop_si=0;iop_si<nop;iop_si++)
	    for(int iflav_so=0;iflav_so<nflavs;iflav_so++)
		  for(int iflav_si=0;iflav_si<nflavs;iflav_si++)
	      {
		int spin_so=meas_pars.mesons[iop_so].first;
		int taste_so=meas_pars.mesons[iop_so].second;
		int spin_si=meas_pars.mesons[iop_si].first;
		int taste_si=meas_pars.mesons[iop_si].second;
		master_fprintf(file," # conf %d ;"
			       " iop_so %d , spin_so %d , taste_so %d ;"
			       " iop_si %d , spin_si %d , taste_si %d ;"
			       " flv_so = %d , m_so = %lg ; flv_si = %d , m_si = %lg\n",
			       iconf,iop_so,spin_so,taste_so,iop_si,spin_si,taste_si,iflav_so,tp.quarks[iflav_so].mass, iflav_si,tp.quarks[iflav_si].mass);
		for(int d_dir=0;d_dir<glbSize[dir];d_dir++)
		  {
		    int ic=icombo(iflav_so,iflav_si,iop_so,iop_si,d_dir,dir);
		    master_fprintf(file,"%d %+16.16lg %+16.16lg\n",d_dir,corr[ic][RE]*norm,corr[ic][IM]*norm);
		  }
		master_fprintf(file,"\n");
	      }
      }
    
    close_file(file);
  }
}
