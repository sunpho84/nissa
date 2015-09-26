#ifndef _NUCLEON_HPP
#define _NUCLEON_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  //parameters to compute time nucleon correlator
  struct nucleon_corr_meas_pars_t
  {
    int flag;
    char path[1024];
    double residue;
    int nhits;
  };
  
  void measure_nucleon_corr(quad_su3 **conf,theory_pars_t theory_pars,nucleon_corr_meas_pars_t meas_pars,int iconf,int conf_created);
}

#endif
