#ifndef _MOMENTA_EVOLVE_HPP
#define _MOMENTA_EVOLVE_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{  
  void accelerate_lx_momenta(quad_su3 *M,quad_su3 *conf,double kappa,int niter,double residue,quad_su3 *H);
  void evolve_lx_momenta_with_force(quad_su3 *H,quad_su3 *F,double dt);
  void evolve_lx_conf_with_momenta(quad_su3 *lx_conf,quad_su3 *H,double dt);
  void evolve_lx_conf_with_accelerated_momenta(quad_su3 *lx_conf,quad_su3 *H,double kappa,int niter,double residue,double dt);
}

#endif
