#ifndef _QUARK_FORCE_HPP
#define _QUARK_FORCE_HPP

#include "multipseudo_rhmc_step.hpp"

namespace nissa
{
  void compute_quark_force_finish_computation(quad_su3 **F,quad_su3 **conf);
  void compute_quark_force_no_stout_remapping(quad_su3 **F,quad_su3 **conf,pseudofermion_t *pf,theory_pars_t *tp,std::vector<rat_approx_t> *appr,std::vector<int> *npfs,double residue);
  void compute_quark_force(quad_su3 **F,quad_su3 **conf,pseudofermion_t *pf,theory_pars_t *physics,std::vector<rat_approx_t> *appr,std::vector<int> *npfs,double residue);

}

#endif
