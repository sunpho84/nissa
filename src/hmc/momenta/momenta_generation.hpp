#ifndef _MOMENTA_GENERATION_HPP
#define _MOMENTA_GENERATION_HPP

#include "new_types/rat_approx.hpp"

namespace nissa
{
  void generate_MFACC_momenta(su3 **pi,int naux_fields,quad_su3 *conf,rat_approx_t *rat_exp_H, double kappa,double residue);
  void generate_hmc_momenta_with_FACC(quad_su3 *H,quad_su3 *conf,rat_approx_t *rat_exp_H,double kappa,double residue);
  void generate_hmc_momenta(quad_su3 **H);
  void generate_hmc_momenta(quad_su3 *H);
}

#endif
