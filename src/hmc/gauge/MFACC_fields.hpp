#ifndef _MFACC_FIELDS_HPP
#define _MFACC_FIELDS_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void generate_MFACC_field(su3 *phi);
  double MFACC_fields_action(su3 **phi,int naux_fields);
  void evolve_MFACC_fields(su3 **phi,int naux_fields,quad_su3 *conf,double kappa,su3 **pi,double dt);
  void evolve_MFACC_momenta(su3 **pi,su3 **phi,int naux_fields,double dt);
  void summ_the_MFACC_momenta_QCD_force(quad_su3 *F,quad_su3 *conf,double kappa,su3 **pi,int naux_fields);
  void summ_the_MFACC_QCD_momenta_QCD_force(quad_su3 *F,quad_su3 *conf,double kappa,int niter,double residue,quad_su3 *H);
}

#endif
