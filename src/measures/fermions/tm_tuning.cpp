#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "tm_tuning.hpp"

namespace nissa
{
  namespace
  {
    const int ncorr_kind=4;
  }
  
  //compute correlation functions for twisted clover, needed to fix tuning
  THREADABLE_FUNCTION_4ARG(tm_tuning, complex*,corr, quad_su3**,conf, quad_u1**,u1b, tm_tuning_meas_pars_t*,meas_pars)
  {
    GET_THREAD_ID();
  }
  THREADABLE_FUNCTION_END
  
  //compute and print
  void measure_tm_tuning(quad_su3 **conf,theory_pars_t &tp,tm_tuning_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    int ncopies=meas_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
        master_fprintf(file,"%d",iconf);
        
        //measure magnetization for each quark
        for(int iflav=0;iflav<tp.nflavs();iflav++)
          {
	    if(tp.quarks[iflav].discretiz!=ferm_discretiz::ROOT_TM_CLOV) crash("not defined for non-Wilson quarks");

	    complex corr[glb_size[0]*ncorr_kind];
	    
            //loop over hits
            int nhits=meas_pars.nhits;
            for(int hit=0;hit<nhits;hit++)
              {
                verbosity_lv2_master_printf("Evaluating tm tuning for flavor %d/%d, ncopies %d/%d nhits %d/%d\n",
                                            iflav+1,tp.nflavs(),icopy+1,ncopies,hit+1,nhits);
            
                //compute and summ
                complex corr_hit[glb_size[0]*ncorr_kind];
                tm_tuning(corr_hit,conf,tp.backfield[iflav],&meas_pars);
                
                //normalize
                for(int i=0;i<ncorr_kind;i++) complex_summ_the_prod_double(corr[i],corr_hit[i],1.0/nhits);
              }
            
            //output
	    for(int t=0;t<glb_size[0];t++)
	      {
		master_fprintf(file,"%d",t);
		for(int ic=0;ic<ncorr_kind;ic++)
		  master_fprintf(file,"\t%+016.16lg \t%+016.16",corr[ic+ncorr_kind*t][RE],corr[ic+ncorr_kind*t][IM]);
		master_fprintf(file,"\n");
	      }
          }
        
        master_fprintf(file,"\n");
      }
    
    close_file(file);
  }
  
  //nucleon correlators
  std::string tm_tuning_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMesonCorrs\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
