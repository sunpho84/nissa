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
  
  //measure the quark number and its derivative w.r.t mu
  void measure_chir_zumba(quad_su3 **conf,theory_pars_t &theory_pars,chir_zumba_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    int nflavs=theory_pars.nflavs();
    
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    complex *point_result=nissa_malloc("point_result",loc_vol,complex);
    NEW_FIELD_T(source);
    
    //vectors for calculation
    NEW_FIELD_T(M);
    NEW_FIELD_T(dM_M);
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	//print conf id
	master_fprintf(file,"%d\t",iconf);
	
	//loop over flavor
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
	    //vectors for output
	    NEW_TRACE_RES_VEC(Tr_M_dM,2);
	    
	    //loop over hits
	    for(int ihit=0;ihit<meas_pars.nhits;ihit++)
	      {
		//fill the source
		fill_source(source);
		
		//compute dM[iflav]*M[jflav]^-1
		MINV(M,iflav,source);
		
		for(int jflav=0;jflav<nflavs;jflav++)
		  {
		    DMDMU(dM_M,iflav,1,M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM[jflav],source,dM_M);
		  }
	      }
	  }
	
	master_fprintf(file,"\n");
      }
    
    DELETE_FIELD_T(M);
  }
  
  //print
  std::string chir_zumba_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasZumba\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }

}
