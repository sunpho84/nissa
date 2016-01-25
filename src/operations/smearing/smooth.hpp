#ifndef _SMOOTH_HPP
#define _SMOOTH_HPP

#include "APE.hpp"
#include "cooling.hpp"
#include "HYP.hpp"
#include "stout.hpp"
#include "Wflow.hpp"

namespace nissa
{
  //parameters to smooth a configuration
  struct smooth_pars_t
  {
    enum method_t{UNSPEC_SMOOTH_METHOD,COOLING,STOUT,WFLOW,APE,HYP};
    
    //basic
    method_t method;
    double meas_each;
    method_t def_method(){return COOLING;}
    double def_meas_each(){return 1;}
    
    //pars
    cool_pars_t cool_pars;
    stout_pars_t stout_pars;
    Wflow_pars_t Wflow_pars;
    ape_pars_t ape_pars;
    hyp_pars_t hyp_pars;
    
    std::string get_method_name()
    {
      switch(method)
	{
	case COOLING:return "Cooling";break;
	case STOUT:return "Stout";break;
	case WFLOW:return "WFlow";break;
	case APE:return "Ape";break;
	case HYP:return "Hyp";break;
	case UNSPEC_SMOOTH_METHOD:return "Unspec";break;
	}
    }
    
    int master_fprintf(FILE *fout,bool full=false);
    int is_nonstandard()
    {
      return
	method!=def_method()||
	meas_each!=def_meas_each();
    }
    
    smooth_pars_t() :
      method(def_method()),
      meas_each(def_meas_each()) {}
  };
}

#endif
