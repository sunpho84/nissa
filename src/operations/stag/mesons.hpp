#ifndef _MESONS_HPP
#define _MESONS_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  struct stag_meson_corr_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int ncopies;
    int nhits;
    
    int def_flag(){return 0;}
    std::string def_path(){return "stag_mesons";}
    double def_residue(){return 1e-12;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    std::vector<std::pair<int,int> > mesons;
    
    int master_fprintf(FILE *fout,bool full=false);
    
    stag_meson_corr_meas_pars_t() : flag(def_flag()),path(def_path()),residue(def_residue()),ncopies(def_ncopies()),nhits(def_nhits()) {}
  };
  
  void measure_staggered_meson_corr(quad_su3 **conf,theory_pars_t &tp,stag_meson_corr_meas_pars_t &pars,int iconf,int conf_created);
}

#endif
