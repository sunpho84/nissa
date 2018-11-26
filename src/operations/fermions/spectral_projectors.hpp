#ifndef _SPECTRAL_PROJECTORS_HPP
#define _SPECTRAL_PROJECTORS_HPP

#include "fermionic_meas.hpp"
#include "hmc/gauge/topological_action.hpp"
#include "hmc/theory_pars.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //parameters to measure topology properties
  struct spectr_proj_meas_pars_t : base_fermionic_meas_t
  {
    int neigs;
    double eig_precision;
    std::string def_path(){return "pettirosso";}
    int def_neigs(){return 5;}
    double def_eig_precision(){return 1e-5;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	path!=def_path() or
	neigs!=def_neigs() or
	eig_precision!=def_eig_precision();
    }
    
    spectr_proj_meas_pars_t() :
      base_fermionic_meas_t(), neigs(def_neigs()), eig_precision(def_eig_precision())
    {path=def_path();}
    virtual ~spectr_proj_meas_pars_t(){}
  };
  
  void measure_spectral_proj(quad_su3 **conf,theory_pars_t &theory_pars,spectr_proj_meas_pars_t &pars, int iconf,bool conf_created);
}

#endif
