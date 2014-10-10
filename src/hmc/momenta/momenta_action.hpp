#ifndef _MOMENTA_ACTION_H
#define _MOMENTA_ACTION_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  double momenta_action(quad_su3 **H);
  double momenta_action(quad_su3 *H);
  double B_momenta_action(double *H_B);
  void MFACC_momenta_action(double *out,su3 **pi,quad_su3 *conf,double kappa);
  inline double MFACC_momenta_action(su3 **pi,quad_su3 *conf,double kappa)
  {
    double out;
    MFACC_momenta_action(&out,pi,conf,kappa);
    return out;
  }
}

#endif
