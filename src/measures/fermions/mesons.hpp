#ifndef _MESONS_HPP
#define _MESONS_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct meson_corr_meas_pars_t :
    base_fermionic_meas_t
  {
    std::vector<std::pair<int,int>> mesons;
    
    std::string def_path() const
    {
      return "meson_corrs";
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false)
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasMesonCorrs\n";
      os<<base_fermionic_meas_t::get_str(full);
      if(mesons.size() or full)
	{
	  os<<" Operators\t=\t{";
	  for(size_t i=0;i<mesons.size();i++)
	    {
	      os<<"("<<mesons[i].first<<","<<mesons[i].second<<")";
	      if(i!=mesons.size()-1) os<<",";
	    }
	  os<<"}\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	mesons.size() or
	path!=def_path();
    }
    
    meson_corr_meas_pars_t() :
      base_fermionic_meas_t()
    {
      path=def_path();
    }
    
    virtual ~meson_corr_meas_pars_t()
    {
    }
  };
  
  void measure_meson_corr(const EoField<quad_su3>& conf,
			  const theory_pars_t& tp,
			  const meson_corr_meas_pars_t& pars,
			  const int& iconf,
			  const int& conf_created);
}

#endif
