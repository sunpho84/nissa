#ifndef _READ_NEW_TYPES_HPP
#define _READ_NEW_TYPES_HPP

#include "operations/smearing/APE.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/smearing/Wflow.hpp"

namespace nissa
{
  void read_Wflow_pars(Wflow_pars_t &pars);
  void read_stout_pars(stout_pars_t &stout_pars);
  void read_ape_pars(ape_pars_t &ape_pars);
}

#endif
