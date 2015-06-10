#ifndef _TIME_PSEUDO_CORR_HPP
#define _TIME_PSEUDO_CORR_HPP

namespace nissa
{
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
  void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created,int dir=0);
}

#endif
