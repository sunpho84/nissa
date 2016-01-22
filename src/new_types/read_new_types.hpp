#ifndef _READ_NEW_TYPES_HPP
#define _READ_NEW_TYPES_HPP

#include "operations/stag/putpourri.hpp"
#include "operations/stag/magnetization.hpp"
#include "operations/stag/mesons.hpp"
#include "operations/stag/nucleon.hpp"
#include "operations/stag/rendens.hpp"
#include "operations/stag/spinpol.hpp"
#include "operations/su3_paths/topological_charge.hpp"

namespace nissa
{
  void read_stout_pars(stout_pars_t &stout_pars);
  void read_ape_pars(ape_pars_t &ape_pars);
}

#endif
