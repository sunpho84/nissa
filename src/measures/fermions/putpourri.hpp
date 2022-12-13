#ifndef _PUTPOURRI_HPP
#define _PUTPOURRI_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct fermionic_putpourri_meas_pars_t : base_fermionic_meas_t
  {
    int compute_susc;
    
    std::string def_path() const
    {
      return "lavanda";
    }
    
    int def_compute_susc() const
    {
      return 0;
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasPutpourri\n";
      if(is_nonstandard() or full)
	os<<base_fermionic_meas_t::get_str(full);
      if(compute_susc!=def_compute_susc() or full)
	os<<" ComputeSusc\t=\t"<<compute_susc<<"\n";
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	compute_susc!=def_compute_susc() or
	path!=def_path();
    }
    
    fermionic_putpourri_meas_pars_t() :
      base_fermionic_meas_t(),
      compute_susc(def_compute_susc())
    {
      path=def_path();
    }
    
    virtual ~fermionic_putpourri_meas_pars_t()
    {
    }
  };
  
  void measure_fermionic_putpourri(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,fermionic_putpourri_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
