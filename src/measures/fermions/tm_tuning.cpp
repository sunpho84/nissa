#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_EXTERNAL_SOLVER
# include "base/export_conf_to_external_solver.hpp"
#endif
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_mix.hpp"
#include "hmc/quark_pars.hpp"
#include "routines/mpi_routines.hpp"

#include "tm_corr_op.hpp"
#include "tm_tuning.hpp"

namespace nissa
{
  /// Compute and print tm tuning
  void measure_tm_tuning(eo_ptr<quad_su3> ext_conf,theory_pars_t &tp,tm_tuning_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    spincolor *eta=nissa_malloc("eta",locVol+bord_vol,spincolor);
    spincolor *phi=nissa_malloc("phi",locVol+bord_vol,spincolor);
    spincolor *phi_r=nissa_malloc("phi_r",locVol+bord_vol,spincolor);
    spincolor *phi_ins_S=nissa_malloc("phi_ins_S",locVol+bord_vol,spincolor);
    spincolor *phi_ins_P=nissa_malloc("phi_ins_P",locVol+bord_vol,spincolor);
    
    /// Store the contractions
    const int ncorr_kind=8;
    complex* contr=nissa_malloc("contr",ncorr_kind*glbSize[0],complex);
    vector_reset(contr);
    
    /// Operations to compute correlators with twisted mass
    tm_corr_op op(ext_conf,meas_pars.residue,tp);
    
    int ncopies=meas_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      for(int iflav=0;iflav<tp.nflavs();iflav++)
	{
	  master_fprintf(file," # conf %d ; flv = %d , m = %lg , k = %lg , cSW = %lg\n",
			 iconf,iflav,tp.quarks[iflav].mass,tp.quarks[iflav].kappa,tp.quarks[iflav].cSW);
	  
	  verbosity_lv2_master_printf("Evaluating tm tuning for flavor %d/%d, ncopies %d/%d\n",
					  iflav+1,tp.nflavs(),icopy+1,ncopies);
	  
	  //Source time
	  coords_t source_coord=generate_random_coord();
	  
	  const int& nhits=meas_pars.nhits;
	  for(int hit=0;hit<nhits;hit++)
	    {
	      //Source time
	      generate_undiluted_source(eta,meas_pars.rnd_type,source_coord[0]);
	      
	      if(hit%2==1)
		{
		  op.inv(phi_r,eta,iflav,1);
#ifdef USE_EXTERNAL_SOLVER
		  export_conf::export_bypass=export_conf::AVOID_EXPORT;
#endif
		}
	      op.inv(phi,eta,iflav,0);
#ifdef USE_EXTERNAL_SOLVER
	      export_conf::export_bypass=export_conf::AVOID_EXPORT;
#endif
	      op.ins(phi_ins_P,5,phi);
	      op.inv(phi_ins_P,phi_ins_P,iflav,0);
	      op.inv(phi_ins_S,phi,iflav,0);
	      if(hit%2==0)
		op.inv(phi_r,eta,iflav,1);
	      
	      auto c=[&](spincolor* oth,int ig,const int icontr)
	      {
		complex temp_contr[glbSize[0]];
		tm_corr_op::undiluted_meson_contr(temp_contr,phi,oth,ig,source_coord[0]);
		for(int t=0;t<glbSize[0];t++)
		  complex_summassign(contr[t+glbSize[0]*icontr],temp_contr[t]);
	      };
	      
	      c(phi,5,0);
	      c(phi_r,5,1);
	      c(phi_r,0,2);
	      c(phi,4,3);
	      c(phi_ins_S,5,4);
	      c(phi_ins_S,4,5);
	      c(phi_ins_P,5,6);
	      c(phi_ins_P,4,7);
	    }
#ifdef USE_EXTERNAL_SOLVER
	  export_conf::export_bypass=export_conf::NO_BYPASS;
#endif	  
	  //output
	  for(int t=0;t<glbSize[0];t++)
	    {
	      master_fprintf(file,"%d  ",t);
	      for(int ic=0;ic<ncorr_kind;ic++)
		{
		  complex c;
		  complex_prod_double(c,contr[t+glbSize[0]*ic],1.0/(meas_pars.nhits*glbSpatVol));
		  master_fprintf(file,"\t%+.16lg , %+.16lg",c[RE],c[IM]);
		}
	      master_fprintf(file,"\n");
	    }
	  
	  master_fprintf(file,"\n");
	}
    
    nissa_free(contr);
    
    nissa_free(eta);
    nissa_free(phi);
    nissa_free(phi_r);
    nissa_free(phi_ins_S);
    nissa_free(phi_ins_P);
    
    close_file(file);
  }
  
  //nucleon correlators
  std::string tm_tuning_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasTmTuning\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
