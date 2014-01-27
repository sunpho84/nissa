#ifndef _STAG_H
#define _STAG_H

namespace nissa
{
  void get_propagator(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
  void measure_chiral_cond(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
  void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created,int dir=0);
}

#endif
