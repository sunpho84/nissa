#ifndef _MAGNETIZATION_HPP
#define _MAGNETIZATION_HPP

namespace nissa
{
  void magnetization(complex *magn,quad_su3 **conf,quark_content_t *quark,color **rnd,color **chi,complex *point_magn,coords *arg,int mu,int nu);
  void magnetization(complex *magn,quad_su3 **conf,int quantization,quad_u1 **u1b,quark_content_t *quark,double residue,color **rnd);
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
}

#endif
