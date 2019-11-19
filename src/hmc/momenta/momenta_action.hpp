#ifndef _MOMENTA_ACTION_HPP
#define _MOMENTA_ACTION_HPP

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  double momenta_action(eo_ptr<quad_su3> H);
  double momenta_action(quad_su3 *H);
  double momenta_action_with_FACC(quad_su3 *conf,double kappa,int niter,double residue,quad_su3 *H);
  void MFACC_momenta_action(double *out,eo_ptr<su3> pi,int naux_fields,quad_su3 *conf,double kappa);
  inline double MFACC_momenta_action(eo_ptr<su3> pi,int naux_fields,quad_su3 *conf,double kappa)
  {
    double out;
    MFACC_momenta_action(&out,pi,naux_fields,conf,kappa);
    return out;
  }
}

#endif
