#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "random/randomGenerate.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "io/input.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"
#include "routines/ios.hpp"

#include "rendens.hpp"
#include "stag.hpp"

namespace nissa
{
#define AT_ORDER(A) if(meas_pars.max_order>=A)
  
  using namespace stag;
  
  //measure the quark number and its derivative w.r.t mu
  void measure_quark_rendens(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,quark_rendens_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    complex *point_result=nissa_malloc("point_result",locVol.nastyConvert(),complex);
    NEW_FIELD_T(source);
    
    //vectors for calculation
    NEW_FIELD_T(M);
    NEW_FIELD_T(dM_M);
    NEW_FIELD_T(d2M_M);
    NEW_FIELD_T(d3M_M);
    NEW_FIELD_T(M_dM_M);
    NEW_FIELD_T(dM_M_dM_M);
    NEW_FIELD_T(d2M_M_dM_M);
    NEW_FIELD_T(M_dM_M_dM_M);
    NEW_FIELD_T(dM_M_dM_M_dM_M);
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	//print conf id
	master_fprintf(file,"%d\t",iconf);
	
	//loop over flavor
	for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
	  {
	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
	    //vectors for output
	    NEW_TRACE_RES(Tr_M_dM);
	    NEW_TRACE_RES(Tr_M_d2M);
	    NEW_TRACE_RES(Tr_M_dM_M_dM);
	    NEW_TRACE_RES(Tr_M_d3M);
	    NEW_TRACE_RES(Tr_M_dM_M_d2M);
	    NEW_TRACE_RES(Tr_M_dM_M_dM_M_dM);
	    
	    //loop over hits
	    for(int ihit=0;ihit<meas_pars.nhits;ihit++)
	      {
		//fill the source
		fill_source(source,-1,meas_pars.rnd_type);
		
		//compute dM*M^-1
		AT_ORDER(1)
		  {
		    MINV(M,iflav,source);
		    DMDMU(dM_M,iflav,1,M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM,source,dM_M);
		  }
		
		//compute d2M*M^-1
		AT_ORDER(2)
		  {
		    DMDMU(d2M_M,iflav,2,M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_d2M,source,d2M_M);
		  }
		
		//compute dM*M^-1*dM*M^-1
		AT_ORDER(2)
		  {
		    MINV(M_dM_M,iflav,dM_M);
		    DMDMU(dM_M_dM_M,iflav,1,M_dM_M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM_M_dM,source,dM_M_dM_M);
		  }
		
		//compute d3M*M^-1
		AT_ORDER(3)
		  {
		    DMDMU(d3M_M,iflav,3,M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_d3M,source,d3M_M);
		  }
		
		//compute d2M*M^-1*dM*M^-1
		AT_ORDER(3)
		  {
		    DMDMU(d2M_M_dM_M,iflav,2,M_dM_M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM_M_d2M,source,d2M_M_dM_M);
		  }
		
		//compute dM_M_dM_M_dM_M
		AT_ORDER(3)
		  {
		    MINV(M_dM_M_dM_M,iflav,dM_M_dM_M);
		    DMDMU(dM_M_dM_M_dM_M,iflav,1,M_dM_M_dM_M);
		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM_M_dM_M_dM,source,dM_M_dM_M_dM_M);
		  }
	      }
	  }
	
	master_fprintf(file,"\n");
      }
    
    DELETE_FIELD_T(M);
    DELETE_FIELD_T(d2M_M);
    DELETE_FIELD_T(dM_M);
    DELETE_FIELD_T(d3M_M);
    DELETE_FIELD_T(M_dM_M);
    DELETE_FIELD_T(dM_M_dM_M);
    DELETE_FIELD_T(d2M_M_dM_M);
    DELETE_FIELD_T(M_dM_M_dM_M);
    DELETE_FIELD_T(dM_M_dM_M_dM_M);
    
    //close and deallocate
    close_file(file);
    nissa_free(point_result);
    DELETE_FIELD_T(source);
  }
  
  //print
  std::string quark_rendens_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasRendens\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(max_order!=def_max_order() or full) os<<" MaxOrder\t=\t"<<max_order<<"\n";
    
    return os.str();
  }
}
