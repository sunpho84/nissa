#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"

#include "stag.hpp"
#include "zumba.hpp"

namespace nissa
{
  using namespace stag;
  
  // measure the chiral condensate and its derivative w.r.t mu
  void measure_chir_zumba(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,chir_zumba_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    crash("reimplement");
    // int nflavs=theory_pars.nflavs();
    
    // //open the file, allocate point result and source
    // FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    // complex *point_result=nissa_malloc("point_result",locVol,complex);
    // NEW_FIELD_T(source);
    
    // //vectors for calculation
    // NEW_FIELD_T(M);           // M^-1
    // NEW_FIELD_T(dM_M);        // M' M^-1
    // NEW_FIELD_T(d2M_M);       // M'' M^-1
    // NEW_FIELD_T(M_M);         // M^-2
    // NEW_FIELD_T(dM_M_M);      // M' M^-2
    // NEW_FIELD_T(d2M_M_M);     // M'' M^-2
    // NEW_FIELD_T(dM_M_dM_M);   // (M' M^-1)^2
    // NEW_FIELD_T(M_dM_M_dM_M); // M^-1 (M' M^-1)^2
    // NEW_FIELD_T(TMP);         // parking variable
    
    // for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
    //   {
    // 	//print conf id and copy id
    // 	master_fprintf(file,"%d\t%d\t",iconf,icopy);
	
    // 	//loop over flavors
    // 	for(int iflav=0;iflav<nflavs;iflav++)
    // 	  {
    // 	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
    // 	    //vectors for output
    // 	    NEW_TRACE_RES(Tr_M);
    // 	    NEW_TRACE_RES(Tr_dM_M);
    // 	    NEW_TRACE_RES(Tr_d2M_M);
    // 	    NEW_TRACE_RES(Tr_M_M);
    // 	    NEW_TRACE_RES(Tr_dM_M_M);
    // 	    NEW_TRACE_RES(Tr_d2M_M_M);
    // 	    NEW_TRACE_RES(Tr_dM_M_dM_M);
    // 	    NEW_TRACE_RES(Tr_M_dM_M_dM_M);
	    
    // 	    //loop over hits
    // 	    for(int ihit=0;ihit<meas_pars.nhits;ihit++)
    // 	      {
    // 		//fill the source
    // 		fill_source(source,-1,meas_pars.rnd_type);
		
    // 		//compute M^-1, M' M^-1, M'' M^-1
    // 		MINV(M,iflav,source);
    // 		DMDMU(dM_M,iflav,1,M);
    // 		DMDMU(d2M_M,iflav,2,M);
		
    // 		//compute M^-2, M' M^-2, M'' M^-2
    // 		MINV(M_M,iflav,M);
    // 		DMDMU(dM_M_M,iflav,1,M_M);
    // 		DMDMU(d2M_M_M,iflav,2,M_M);
		
    // 		//compute (M' M^-1)^2
    // 		DMDMU(dM_M_dM_M,iflav,1,M); // M' M^-1
    // 		MINV(TMP,iflav,dM_M_dM_M); // M^-1 M' M^-1
    // 		DMDMU(dM_M_dM_M,iflav,1,TMP);
		
    // 		//compute M^-1 (M' M^-1)^2
    // 		MINV(M_dM_M_dM_M,iflav,dM_M_dM_M);
		
    // 		//print traces
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M,source,M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_dM_M,source,dM_M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_d2M_M,source,d2M_M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_M,source,M_M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_dM_M_M,source,dM_M_M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_d2M_M_M,source,d2M_M_M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_dM_M_dM_M,source,dM_M_dM_M);
    // 		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM_M_dM_M,source,M_dM_M_dM_M);
    // 	      }
    // 	  }
	
    // 	master_fprintf(file,"\n");
    //   }
    
    // //deallocate and close file
    // DELETE_FIELD_T(M);
    // DELETE_FIELD_T(dM_M);
    // DELETE_FIELD_T(d2M_M);
    // DELETE_FIELD_T(M_M);
    // DELETE_FIELD_T(dM_M_M);
    // DELETE_FIELD_T(d2M_M_M);
    // DELETE_FIELD_T(dM_M_dM_M);
    // DELETE_FIELD_T(M_dM_M_dM_M);
    // DELETE_FIELD_T(TMP);
    
    // close_file(file);
    // nissa_free(point_result);
    // DELETE_FIELD_T(source);
  }
}
