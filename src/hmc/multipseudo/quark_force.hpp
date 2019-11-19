#ifndef _QUARK_FORCE_HPP
#define _QUARK_FORCE_HPP

#include "multipseudo_rhmc_step.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void compute_quark_force_finish_computation(eo_ptr<quad_su3> F,eo_ptr<quad_su3> conf);
  void compute_quark_force_no_stout_remapping(eo_ptr<quad_su3> F,eo_ptr<quad_su3> conf,std::vector<std::vector<pseudofermion_t> > *pf,theory_pars_t *tp,std::vector<rat_approx_t> *appr,double residue);
  void compute_quark_force(eo_ptr<quad_su3> F,eo_ptr<quad_su3> conf,std::vector<std::vector<pseudofermion_t> > *pf,theory_pars_t *physics,std::vector<rat_approx_t> *appr,double residue);

}

#endif
