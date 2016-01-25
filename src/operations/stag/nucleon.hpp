#ifndef _NUCLEON_HPP
#define _NUCLEON_HPP

#include "stag.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  //parameters to compute time nucleon correlator
  struct nucleon_corr_meas_pars_t : base_fermionic_meas_t
  {
    std::string def_path(){return "nucleon_corr";}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	path!=def_path();
    }
    
    nucleon_corr_meas_pars_t() :
      base_fermionic_meas_t()
    {path=def_path();}
    virtual ~nucleon_corr_meas_pars_t(){}
  };
  
  void measure_nucleon_corr(quad_su3 **conf,theory_pars_t theory_pars,nucleon_corr_meas_pars_t meas_pars,int iconf,int conf_created);
}

#endif
