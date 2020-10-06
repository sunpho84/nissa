#ifndef _MESONS_HPP
#define _MESONS_HPP

#include "stag.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct meson_corr_meas_pars_t : base_fermionic_meas_t
  {
    std::vector<std::pair<int,int> > mesons;
    
    std::string def_path(){return "meson_corrs";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	mesons.size() or
	path!=def_path();
    }
    
    meson_corr_meas_pars_t() :
      base_fermionic_meas_t()
    {path=def_path();}
    virtual ~meson_corr_meas_pars_t(){}
  };
  
  void measure_meson_corr(quad_su3 **conf,theory_pars_t &tp,meson_corr_meas_pars_t &pars,int iconf,int conf_created);
}

#endif
