#ifndef _MESONS_HPP
#define _MESONS_HPP

#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct meson_corr_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    double residue;
    int ncopies;
    int nhits;
    std::vector<std::pair<int,int> > mesons;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "mesons";}
    double def_residue(){return 1e-12;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	path!=def_path()||
	residue!=def_residue()||
	ncopies!=def_ncopies()||
	nhits!=def_nhits();
    }
    
    meson_corr_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      residue(def_residue()),
      ncopies(def_ncopies()),
      nhits(def_nhits()) {}
  };
  
  void measure_staggered_meson_corr(quad_su3 **conf,theory_pars_t &tp,meson_corr_meas_pars_t &pars,int iconf,int conf_created);
}

#endif
