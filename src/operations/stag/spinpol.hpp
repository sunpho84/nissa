#ifndef _SPINPOL_HPP
#define _SPINPOL_HPP

#include "stag.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  struct spinpol_meas_pars_t : base_fermionic_meas_t
  {
    int dir;
    int use_ferm_conf_for_gluons;
    smooth_pars_t smooth_pars;
    
    std::string def_path(){return "pollo";}
    int def_dir(){return 1;}
    int def_use_ferm_conf_for_gluons(){return 0;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	dir!=def_dir()||
	use_ferm_conf_for_gluons!=def_use_ferm_conf_for_gluons()||
	path!=def_path()||
	smooth_pars.is_nonstandard();
    }
    
    spinpol_meas_pars_t() :
      base_fermionic_meas_t(),
      dir(def_dir()),
      use_ferm_conf_for_gluons(def_use_ferm_conf_for_gluons())
    {path=def_path();}
  };
  
  void measure_spinpol(quad_su3 **ferm_conf,quad_su3 **glu_conf,theory_pars_t &theory_pars,spinpol_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
