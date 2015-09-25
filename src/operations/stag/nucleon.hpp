#ifndef _NUCLEON_HPP
#define _NUCLEON_HPP

namespace nissa
{
  //parameters to compute time nucleon correlator
  struct nucleon_meas_pars_t
  {
    int flag;
    char path[1024];
    double residue;
    int nhits;
  };
}

#endif
