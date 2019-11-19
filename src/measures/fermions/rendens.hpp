#ifndef _RENDENS_HPP
#define _RENDENS_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct quark_rendens_meas_pars_t : base_fermionic_meas_t
  {
    int max_order;
    
    int def_max_order(){return 2;}
    std::string def_path(){return "rende";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	max_order!=def_max_order()||
	path!=def_path();
    }
    
    quark_rendens_meas_pars_t() :
      base_fermionic_meas_t(),
      max_order(def_max_order())
    {path=def_path();}
    virtual ~quark_rendens_meas_pars_t(){}
  };
  
  void measure_quark_rendens(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,quark_rendens_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
