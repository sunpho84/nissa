#ifndef _ZUMBA_HPP
#define _ZUMBA_HPP

#include "hmc/theory_pars.hpp"

#include "fermionic_meas.hpp"

namespace nissa
{
  struct chir_zumba_meas_pars_t :
    base_fermionic_meas_t
  {
    int max_order;
    
    std::string def_path() const
    {
      return "zambia";
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasZumba\n";
      os<<base_fermionic_meas_t::get_str(full);
      
      return os.str();
    }
    
    bool is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path();
    }
    
    chir_zumba_meas_pars_t() :
      base_fermionic_meas_t()
    {
      path=def_path();
    }
    
    virtual ~chir_zumba_meas_pars_t()
    {
    }
  };
  
  void measure_chir_zumba(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,chir_zumba_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
