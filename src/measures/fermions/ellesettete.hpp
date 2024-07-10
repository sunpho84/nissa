#ifndef _ELLESETTETE_HPP
#define _ELLESETTETE_HPP

#include "hmc/theory_pars.hpp"

#include "fermionic_meas.hpp"

namespace nissa
{
  struct ellesettete_meas_pars_t : base_fermionic_meas_t
  {
    int max_order;
    int method; // true for numerical, false for analytical
    double epsilon; // precision for numerical method
    std::string def_path(){return "ellesettete";}
    int def_method(){return 0;}
    double def_epsilon(){return 1e-6;}
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path() or
  method!=def_method() or
  epsilon!=def_epsilon();
    }
    
    ellesettete_meas_pars_t() :
      base_fermionic_meas_t(),
      method(def_method()),
      epsilon(def_epsilon())
    {
    path=def_path();
    }
    virtual ~ellesettete_meas_pars_t(){}
  };
  
  void measure_ellesettete(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,ellesettete_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
