#ifndef _MINMAX_EIGENVALUES_HPP
#define _MINMAX_EIGENVALUES_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct minmax_eigenvalues_meas_pars_t : base_fermionic_meas_t
  {
    std::string def_path(){return "plover";}
    
    int neigs;
    int wspace_size;
    int min_max;
    int def_neigs(){return 5;}
    int def_min_max(){return 0;}
    int def_wspace_size(){return 100;}
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path() or
        neigs!=def_neigs() or
        wspace_size!=def_wspace_size() or
        min_max!=def_min_max();
    }
    
    minmax_eigenvalues_meas_pars_t() :
      base_fermionic_meas_t(),
      neigs(def_neigs()),
      wspace_size(def_wspace_size()),
      min_max(def_min_max())
    {
      path=def_path();
    }
    virtual ~minmax_eigenvalues_meas_pars_t(){}
  };
  void measure_minmax_eigenvalues(quad_su3 **conf,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
