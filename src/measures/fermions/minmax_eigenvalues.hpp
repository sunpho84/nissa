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
    double eig_precision;
    int def_neigs(){return 5;}
    double def_eig_precision(){return 1e-5;}
    int min_max;
    int def_min_max(){return 0;}
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	path!=def_path() or
        neigs!=def_neigs() or
        eig_precision!=def_eig_precision() or
        min_max!=def_min_max();
    }
    
    minmax_eigenvalues_meas_pars_t() :
    base_fermionic_meas_t(),
    neigs(def_neigs()),
    eig_precision(def_eig_precision()),
    min_max(def_min_max())
    {
      path=def_path();
    }
    virtual ~minmax_eigenvalues_meas_pars_t(){}
  };
  void measure_minmax_eigenvalues(quad_su3 **conf,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
