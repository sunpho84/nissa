#ifndef _SPINPOL_HPP
#define _SPINPOL_HPP

#include "operations/smearing/smooth.hpp"

namespace nissa
{
  struct spinpol_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    double residue;
    int dir;
    int nhits;
    int use_ferm_conf_for_gluons;
    smooth_pars_t smooth_pars;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "pollo";}
    double def_residue(){return 1e-12;}
    int def_dir(){return 1;}
    int def_nhits(){return 1;}
    int def_use_ferm_conf_for_gluons(){return 0;}
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	path!=def_path()||
	residue!=def_residue()||
	dir!=def_dir()||
	nhits!=def_nhits()||
	use_ferm_conf_for_gluons!=def_use_ferm_conf_for_gluons()||
	smooth_pars.is_nonstandard();
    }
    
    spinpol_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      residue(def_residue()),
      dir(def_dir()),
      nhits(def_nhits()),
      use_ferm_conf_for_gluons(def_use_ferm_conf_for_gluons()) {}
  };
  
  void measure_spinpol(quad_su3 **ferm_conf,quad_su3 **glu_conf,theory_pars_t &theory_pars,spinpol_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
