#ifndef _RENDENS_HPP
#define _RENDENS_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct quark_rendens_meas_pars_t :
    base_fermionic_meas_t
  {
    int max_order;
    
    int def_max_order() const
    {
      return 2;
    }
    
    std::string def_path() const
    {
      return "rende";
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasRendens\n";
      os<<base_fermionic_meas_t::get_str(full);
      if(max_order!=def_max_order() or full) os<<" MaxOrder\t=\t"<<max_order<<"\n";
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	max_order!=def_max_order() or
	path!=def_path();
    }
    
    quark_rendens_meas_pars_t() :
      base_fermionic_meas_t(),
      max_order(def_max_order())
    {
      path=def_path();
    }
    
    virtual ~quark_rendens_meas_pars_t()
    {
    }
  };
  
  void measure_quark_rendens(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,quark_rendens_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
