#ifndef _STAG_HPP
#define _STAG_HPP

namespace nissa
{
  void magnetization(complex *magn,quad_su3 **conf,quark_content_t *quark,color **rnd,color **chi,complex *point_magn,coords *arg,int mu,int nu);
  void magnetization(complex *magn,quad_su3 **conf,int quantization,quad_u1 **u1b,quark_content_t *quark,double residue,color **rnd);
  void get_propagator(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
  void measure_fermionic_putpourri(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
  void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created,int dir=0);
}

#endif
