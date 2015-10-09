#ifndef _NUCLEON_HPP
#define _NUCLEON_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  //parameters to compute time nucleon correlator
  struct nucleon_corr_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int nhits;
    
    void master_fprintf(FILE *fout);
    
    nucleon_corr_meas_pars_t() : flag(0),path("nucleon_corr"),residue(1e-12),nhits(1) {}
  };
  
  void measure_nucleon_corr(quad_su3 **conf,theory_pars_t theory_pars,nucleon_corr_meas_pars_t meas_pars,int iconf,int conf_created);
}

#endif
