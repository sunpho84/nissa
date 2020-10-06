#ifndef _ZUMBA_HPP
#define _ZUMBA_HPP

#include "hmc/theory_pars.hpp"

#include "stag.hpp"

namespace nissa
{
  struct chir_zumba_meas_pars_t : base_fermionic_meas_t
  {
    int max_order;
    
    std::string def_path(){return "zambia";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path();
    }
    
    chir_zumba_meas_pars_t() :
      base_fermionic_meas_t()
    {path=def_path();}
    virtual ~chir_zumba_meas_pars_t(){}
  };
  
  void measure_chir_zumba(quad_su3 **conf,theory_pars_t &theory_pars,chir_zumba_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
