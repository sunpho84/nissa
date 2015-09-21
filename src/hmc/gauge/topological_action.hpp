#ifndef _TOPOLOGICAL_ACTION_HPP
#define _TOPOLOGICAL_ACTION_HPP

namespace nissa
{
  double topodynamical_potential(double Q,topotential_pars_t &pars);
  void save_topodynamical_potential(topotential_pars_t &pars);
  void load_topodynamical_potential(topotential_pars_t &pars);
  double topotential_action(quad_su3 **ext_conf,topotential_pars_t &pars);
}

#endif
