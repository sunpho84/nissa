#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "random/randomGenerate.hpp"
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
    
    spincolor *eta=nissa_malloc("eta",locVolWithBord.nastyConvert(),spincolor);
    spincolor *phi=nissa_malloc("phi",locVolWithBord.nastyConvert(),spincolor);
    spincolor *phi_ins_S=nissa_malloc("phi_ins_S",locVolWithBord.nastyConvert(),spincolor);
    spincolor *phi_ins_P=nissa_malloc("phi_ins_P",locVolWithBord.nastyConvert(),spincolor);
    
    /// Store the contractions
    const int ncorr_kind=6;
    complex* contr=nissa_malloc("contr",ncorr_kind*glbTimeSize.nastyConvert(),complex);
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
	  GlbCoords source_coord;
	  generate_random_coord(source_coord);
	  
	  const int& nhits=meas_pars.nhits;
	  for(int hit=0;hit<nhits;hit++)
	    {
	      //Source time
	      generate_undiluted_source(eta,meas_pars.rnd_type,source_coord(timeDirection));
	      
	      op.inv(phi,eta,iflav);
	      op.ins(phi_ins_P,5,phi);
	      op.inv(phi_ins_P,phi_ins_P,iflav);
	      op.inv(phi_ins_S,phi,iflav);
	      
	      auto c=[&](spincolor* oth,int ig,const int icontr)
	      {
		complex* temp_contr=new complex[glbTimeSize.nastyConvert()];
		tm_corr_op::undiluted_meson_contr(temp_contr,phi,oth,ig,source_coord(timeDirection));
		FOR_ALL_GLB_TIMES(t)
		  complex_summassign(contr[t.nastyConvert()+glbTimeSize.nastyConvert()*icontr],temp_contr[t.nastyConvert()]);
		
		delete[] temp_contr;
	      };
	      
	      c(phi,5,0);
	      c(phi,4,1);
	      c(phi_ins_S,5,2);
	      c(phi_ins_S,4,3);
	      c(phi_ins_P,5,4);
	      c(phi_ins_P,4,5);
	    }
	  
	  //output
	  FOR_ALL_GLB_TIMES(t)
	    {
	      master_fprintf(file,"%d  ",t);
	      for(int ic=0;ic<ncorr_kind;ic++)
		{
		  complex c;
		  complex_prod_double(c,contr[t.nastyConvert()+glbTimeSize.nastyConvert()*ic],1.0/(meas_pars.nhits*glbSpatVol()));
		  master_fprintf(file,"\t%+.16lg , %+.16lg",c[RE],c[IM]);
		}
	      master_fprintf(file,"\n");
	    }
	  
	  master_fprintf(file,"\n");
	}
    
    nissa_free(contr);
    
    nissa_free(eta);
    nissa_free(phi);
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
