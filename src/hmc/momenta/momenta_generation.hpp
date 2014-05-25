#ifndef _MOMENTA_GENERATION_H
#define _MOMENTA_GENERATION_H

namespace nissa
{
  void generate_hmc_momenta(quad_su3 **H);
  void generate_hmc_B_momenta(double *H_B);
}

#endif
