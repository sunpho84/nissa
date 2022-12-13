#ifndef _TM_TUNING_HPP
#define _TM_TUNING_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct tm_tuning_meas_pars_t :
    base_fermionic_meas_t
  {
    std::string def_path() const
    {
      return "tm_tuning";
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasTmTuning\n";
      os<<base_fermionic_meas_t::get_str(full);
      
      return os.str();
    }
    
    bool is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path();
    }
    
    tm_tuning_meas_pars_t() :
      base_fermionic_meas_t()
    {
      path=def_path();
    }
    
    virtual ~tm_tuning_meas_pars_t()
    {
    }
  };
  
  void measure_tm_tuning(eo_ptr<quad_su3> conf,theory_pars_t &tp,tm_tuning_meas_pars_t &pars,int iconf,int conf_created);
}

#endif
