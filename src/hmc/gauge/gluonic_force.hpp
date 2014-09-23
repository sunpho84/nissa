#ifndef _GLUONIC_FORCE_H
#define _GLUONIC_FORCE_H

namespace nissa
{
  void gluonic_force_finish_computation(quad_su3*F,quad_su3*conf,bool phase_pres);
  void compute_gluonic_force_lx_conf(quad_su3*F,quad_su3*conf,theory_pars_t*physics,bool phase_pres);
}

#endif
