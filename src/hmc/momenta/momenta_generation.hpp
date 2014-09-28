#ifndef _MOMENTA_GENERATION_HPP
#define _MOMENTA_GENERATION_HPP

namespace nissa
{
  void generate_MFACC_momenta(su3 **pi,quad_su3 *conf,double kappa,double residue);
  void generate_hmc_momenta(quad_su3 *H);
  void generate_hmc_momenta(quad_su3 **H);
  void generate_hmc_B_momenta(double *H_B);
}

#endif
