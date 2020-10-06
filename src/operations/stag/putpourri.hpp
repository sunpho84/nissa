#ifndef _PUTPOURRI_HPP
#define _PUTPOURRI_HPP

#include "stag.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct fermionic_putpourri_meas_pars_t : base_fermionic_meas_t
  {
    int compute_susc;
    
    std::string def_path(){return "lavanda";}
    int def_compute_susc(){return 0;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	compute_susc!=def_compute_susc()||
	path!=def_path();
    }
    
    fermionic_putpourri_meas_pars_t() :
      base_fermionic_meas_t(),
      compute_susc(def_compute_susc())
    {path=def_path();}
    virtual ~fermionic_putpourri_meas_pars_t(){}
  };
  
  void measure_fermionic_putpourri(quad_su3 **conf,theory_pars_t &theory_pars,fermionic_putpourri_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
