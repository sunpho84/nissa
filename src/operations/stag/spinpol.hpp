#ifndef _SPINPOL_HPP
#define _SPINPOL_HPP

namespace nissa
{
  void measure_spinpol(quad_su3 **ferm_conf,quad_su3 **glu_conf,theory_pars_t &theory_pars,spinpol_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
